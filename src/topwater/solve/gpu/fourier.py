"""GPU Fourier synthesis with unchanged native lengths, including odd sizes.

Bluestein reduces arbitrary DFT lengths to a radix-2 convolution. No resampling,
padding of physical motion, or changed output clock is introduced; convolution
padding is internal to the exact DFT algorithm.
"""
from importlib import resources
import numpy as np
from scipy.fft import fft

from topwater.solve.fourier import FourierPreparation
from topwater.solve.transfer import Transfer, prepare_host_transfer
from topwater.solve.gpu.memory import ReusableBuffers


class GPUFourierWorkspace(FourierPreparation):
    def __init__(self,evaluator):
        super().__init__()
        self.evaluator=evaluator
        self.cl,self.context,self.queue=evaluator.cl,evaluator.context,evaluator.queue
        code=resources.files(__package__).joinpath('fourier.cl').read_text(encoding='utf-8')
        self.program=self.cl.Program(self.context,code).build()
        for name in ('fft_stage','start_real','start_inverse','convolve','finish_forward',
                     'finish_selected','multiply_spectrum'):
            setattr(self,name,self.cl.Kernel(self.program,name))
        self.n=None
        self.capacity=0
        self.buffers=ReusableBuffers(self.queue)

    def upload(self,values):
        values=np.ascontiguousarray(values)
        return self.cl.Buffer(self.context,self.cl.mem_flags.READ_ONLY|self.cl.mem_flags.COPY_HOST_PTR,
                              hostbuf=values)

    def buffer(self,size):
        return self.cl.Buffer(self.context,self.cl.mem_flags.READ_WRITE,size)

    def plan(self,n):
        if self.n==n:return
        self.n=n
        self.m=1<<max(1,(2*n-2).bit_length())
        self.bits=self.m.bit_length()-1
        j=np.arange(n,dtype=np.int64)
        chirp=np.exp(-1j*np.pi*(j*j%(2*n))/n)
        kernel=np.zeros(self.m,complex)
        kernel[:n]=chirp.conj()
        kernel[self.m-j[1:]]=chirp[1:].conj()
        self.chirp=self.upload(chirp)
        self.convolution_kernel=self.upload(fft(kernel))
        self.twiddle=self.upload(np.exp(-2j*np.pi*np.arange(self.m//2)/self.m))
        self.capacity=0

    def work(self,rows):
        if rows>self.capacity:
            size=rows*self.m*16
            self.buffers.get('work_a',size)
            self.buffers.get('work_b',size)
            self.capacity=rows

    def reserve(self,n,rows,output_count):
        self.plan(n)
        self.work(max(6,rows))
        frequencies=n//2+1
        for name,size in dict(spectrum=rows*frequencies*16, result=rows*output_count*8,
                initial=rows*8, premax=rows*8, selected=output_count*4,
                velocity=6*frequencies*16, velocity_values=6*n*8).items():
            self.buffers.get(name,size)

    def reset_buffers(self):
        self.buffers.clear()
        self.capacity=0
        self.reset_preparation()

    @property
    def nbytes(self):
        return self.buffers.nbytes

    def copy(self,name,values):
        values=np.ascontiguousarray(values)
        target=self.buffers.get(name,values.nbytes)
        self.cl.enqueue_copy(self.queue,target,values).wait()
        return target

    def power_two(self,buffer,rows,inverse):
        for span in (1<<power for power in range(1,self.bits+1)):
            self.fft_stage(self.queue,(rows*self.m//2,),None,buffer,np.int32(self.m),
                           np.int32(rows),np.int32(span),np.int32(inverse),self.twiddle)

    def convolution(self,rows):
        self.power_two(self.buffers['work_a'],rows,False)
        self.convolve(self.queue,(rows*self.m,),None,self.buffers['work_a'],self.convolution_kernel,
            np.int32(self.m),np.int32(rows),np.int32(self.bits),self.buffers['work_b'])
        self.power_two(self.buffers['work_b'],rows,True)

    def prepare_transform(self,values,n):
        self.plan(n)
        self.work(6)
        value_buffer=self.copy('velocity_values',np.ascontiguousarray(values.T))
        self.start_real(self.queue,(6*self.m,),None,value_buffer,np.int32(n),np.int32(self.m),
                        np.int32(6),np.int32(self.bits),self.chirp,self.buffers['work_a'])
        self.convolution(6)
        velocity_buffer=self.buffers.get('velocity',6*len(self.f)*16)
        self.finish_forward(self.queue,(6*len(self.f),),None,self.buffers['work_b'],self.chirp,
            np.int32(n),np.int32(self.m),np.int32(6),velocity_buffer)

    def prepare_transfer(self,transfer,nw):
        if isinstance(transfer,Transfer):
            if transfer.context!=self.context or transfer.shape[0]!=nw or transfer.shape[-1]!=6:
                raise ValueError('GPU transfer belongs to a different clock, device, or mode layout')
            return transfer
        return prepare_host_transfer(transfer,nw)

    def spectrum(self,H,w,values,n,dt):
        self.prepare(w,values,n,dt)
        rows=H.shape[1]
        spectrum_buffer=self.buffers.get('spectrum',rows*len(self.f)*16)
        device=isinstance(H,Transfer)
        coefficients=H.buffer if device else self.upload(np.ascontiguousarray(H))
        self.multiply_spectrum(self.queue,(rows*len(self.f),),None,coefficients,
            np.int32(rows),self.buffers['velocity'],np.int32(len(self.f)),
            np.int32(H.leading_zero if device else 0),np.float64(H.scale if device else 1.),spectrum_buffer)
        return rows

    def inverse_response(self,rows,n,zero,selected):
        self.work(rows)
        self.start_inverse(self.queue,(rows*self.m,),None,self.buffers['spectrum'],np.int32(n),
            np.int32(self.m),np.int32(rows),np.int32(self.bits),self.chirp,self.buffers['work_a'])
        self.convolution(rows)
        result=np.empty((len(selected),rows),float)
        initial=np.empty(rows,float);premax=np.empty(rows,float)
        result_buffer=self.buffers.get('result',result.nbytes)
        initial_buffer=self.buffers.get('initial',initial.nbytes);premax_buffer=self.buffers.get('premax',premax.nbytes)
        selected_buffer=self.copy('selected',np.asarray(selected,dtype=np.int32))
        self.finish_selected(self.queue,(rows*len(selected),),None,self.buffers['work_b'],self.chirp,
            np.int32(n),np.int32(self.m),np.int32(rows),np.int32(zero),
            selected_buffer,np.int32(len(selected)),
            result_buffer,initial_buffer,premax_buffer)
        self.cl.enqueue_copy(self.queue,result,result_buffer).wait()
        self.cl.enqueue_copy(self.queue,initial,initial_buffer).wait()
        self.cl.enqueue_copy(self.queue,premax,premax_buffer).wait()
        return result,initial.tolist(),premax.tolist() if zero else None

    def close(self):
        try:
            self.reset_buffers()
        finally:
            for name,value in list(vars(self).items()):
                if isinstance(value,self.cl.Buffer):
                    value.release();setattr(self,name,None)
            self.program=None
