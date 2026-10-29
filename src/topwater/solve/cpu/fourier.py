"""Reusable host Fourier buffers and CPU numerical operations."""
import numpy as np
from scipy.fft import rfft
from topwater.solve.fourier import FourierPreparation
from topwater.solve.transfer import Transfer, prepare_host_transfer


class CPUFourierWorkspace(FourierPreparation):
    def __init__(self):
        super().__init__()
        self._spectrum = self._field = None

    prepare_transfer = staticmethod(prepare_host_transfer)

    def reserve(self,n,rows,output_count):
        shape=(n//2+1,rows)
        if (self._spectrum is None or self._spectrum.shape[0] != shape[0]
                or self._spectrum.shape[1] < rows or self._field.shape[0] != n):
            self._spectrum=np.empty(shape,complex)
            self._field=np.empty((n,rows),float)
        self.m=n

    @property
    def nbytes(self):
        return sum(value.nbytes for value in (self._spectrum,self._field) if value is not None)

    def reset_buffers(self):
        self._spectrum=self._field=None
        self.reset_preparation()
        self.velocity_fft=None

    close=reset_buffers

    def prepare_transform(self,values,n):
        self.velocity_fft = rfft(values,axis=0)

    def spectrum(self,H,w,values,n,dt):
        self.prepare(w,values,n,dt)
        self.reserve(n,H.shape[1],n)
        spectrum = self._spectrum[:,:H.shape[1]]
        spectrum.fill(0.)
        def gather(indices,mode):
            if not isinstance(H,Transfer):
                return H[indices,:,mode]
            rows=H.buffer[np.maximum(indices-int(H.leading_zero),0),:,mode]
            if H.leading_zero:
                rows[indices==0]=0.
            return rows*H.scale
        indices = np.arange(len(w))
        for mode in range(6):
            spectrum += gather(indices,mode).conj()*self.velocity_fft[:,mode,None]
        return spectrum

    def inverse_response(self,spectrum,n,zero,selected):
        field = self._field[:,:spectrum.shape[1]]
        np.fft.irfft(spectrum,n=n,axis=0,out=field)
        pre_max = np.max(abs(field[:zero]),axis=0).tolist() if zero else None
        return field[zero+selected].copy(),field[zero].tolist(),pre_max
