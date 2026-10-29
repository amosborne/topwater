"""Potential-only OpenCL evaluation for the mean-water-surface workflow.

Supports Capytaine 3.0.0, infinite depth, z=0 targets, low_freq singularities,
scaled_nemoh3 tabulation and float64 only. It does not solve body BEM equations.
PyOpenCL is imported only when GPU reconstruction is selected.
"""
from importlib import resources
from time import perf_counter

import numpy as np
from topwater.solve.gpu.memory import ReusableBuffers, GPUAllocationError
from topwater.solve.transfer import Transfer as DeviceTransfer


class OpenCLPotential:
    allocation_error = GPUAllocationError
    def __init__(self, green_function, source_mesh, *, device_name=None):
        import pyopencl as cl
        from topwater.solve.capytaine import potential_inputs
        self.cl = cl
        inputs = potential_inputs(green_function, source_mesh)
        available = [device for platform in cl.get_platforms() for device in platform.get_devices()
                     if device.type & cl.device_type.GPU and device.double_fp_config]
        candidates = available
        if device_name is not None:
            exact = [device for device in available if device.name.casefold() == device_name.casefold()]
            candidates = exact or [device for device in available if device_name.casefold() in device.name.casefold()]
            if len({device.name.casefold() for device in candidates}) > 1:
                raise RuntimeError('Ambiguous GPU name; matching devices: '+', '.join(dict.fromkeys(device.name for device in candidates)))
        if not candidates:
            names = ', '.join(dict.fromkeys(device.name for device in available)) or 'none'
            raise RuntimeError(f'No double-precision OpenCL GPU matching {device_name!r}; available GPUs: {names}')
        # Prefer dedicated devices; integrated GPUs can report more shared RAM.
        self.device = max(candidates, key=lambda d: (not d.host_unified_memory,
                                                    d.max_compute_units, d.global_mem_size))
        self.context = cl.Context([self.device])
        self.queue = cl.CommandQueue(self.context, properties=cl.command_queue_properties.PROFILING_ENABLE)
        self._closed = False
        try:
            self._initialize(inputs)
        except BaseException as error:
            try:
                self.close()
            except Exception as cleanup_error:
                error.add_note(f'GPU evaluator initialization cleanup also failed: {cleanup_error}')
            raise

    def _initialize(self, inputs):
        cl = self.cl
        self.batch_buffers = ReusableBuffers(self.queue)
        started = perf_counter()
        code = resources.files(__package__).joinpath('potential.cl').read_text(encoding='utf-8')
        self.program = cl.Program(self.context, code).build()
        self.build_seconds = perf_counter()-started
        self.rankine_kernel = cl.Kernel(self.program, 'rankine_matrix')
        self.potential_kernel = cl.Kernel(self.program, 'potential_matrix')
        self.incident_kernel = cl.Kernel(self.program, 'incident_spectrum')
        self.apply_kernel = cl.Kernel(self.program, 'apply_sources')
        self.pi, self.gamma = inputs.pi, inputs.gamma
        self.nf = len(inputs.centers)
        self.nq = inputs.weights.shape[1]
        self.nr, self.nz = inputs.table.shape[:2]
        for name in ('nodes', 'centers', 'normals', 'areas', 'radii', 'quadrature', 'weights', 'rr', 'zz', 'table'):
            setattr(self, name, self.upload(getattr(inputs, name)))
        self._points = None
        self._source_values = None

    def upload(self, array):
        values = np.ascontiguousarray(array, dtype=np.float64)
        return self.cl.Buffer(self.context, self.cl.mem_flags.READ_ONLY | self.cl.mem_flags.COPY_HOST_PTR,
                              hostbuf=values)

    def close(self):
        """Release per-solve resources; safe to call repeatedly."""
        if getattr(self, '_closed', False):
            return
        try:
            if hasattr(self, 'batch_buffers'):
                self.batch_buffers.clear()
            else:
                self.queue.finish()
        finally:
            for name, value in list(vars(self).items()):
                if isinstance(value, self.cl.Buffer):
                    value.release()
                    setattr(self, name, None)
            self.rankine_kernel = self.potential_kernel = self.incident_kernel = self.apply_kernel = None
            self._source_values = None
            self.program = self.queue = self.context = None
            self._closed = True

    def set_points(self, points):
        if self._closed:
            raise RuntimeError('GPU field evaluator is closed')
        values = np.asarray(points, dtype=float)
        if (values.ndim != 2 or values.shape[1] != 3 or not len(values)
                or not np.isfinite(values).all() or np.any(values[:, 2] != 0)):
            raise ValueError('Prototype requires finite nonempty (point,3) z=0 targets')
        if self._points is not None and np.array_equal(values, self._points):
            return
        self._points = None
        points_buffer = self.batch_buffers.get('points', values.nbytes)
        rankine = self.batch_buffers.get('rankine', len(values)*self.nf*8)
        self.cl.enqueue_copy(self.queue, points_buffer, np.ascontiguousarray(values)).wait()
        self.rankine_kernel(self.queue, (len(values)*self.nf,), None,
            points_buffer, np.int32(len(values)), np.int32(self.nf), self.nodes, self.centers,
            self.normals, self.areas, self.radii, self.pi, rankine).wait()
        self._points, self.np = values.copy(), len(values)

    def matrix_budget(self):
        return min(256*1024**2, self.device.max_mem_alloc_size, self.device.global_mem_size//8)

    def reserve_radiation(self, rows, frequencies):
        self._points = self._source_values = None
        budget = self.matrix_budget()
        points = min(rows, max(1, budget//(self.nf*16)))
        for name, size in dict(points=points*3*8, rankine=points*self.nf*8,
                matrix=min(frequencies*points*self.nf*16, budget), transfer=frequencies*rows*6*16,
                sources=frequencies*self.nf*6*16, wavenumbers=frequencies*8).items():
            self.batch_buffers.get(name, size)

    def reset_batch_buffers(self):
        self.batch_buffers.clear()
        self._points = self._source_values = None

    @property
    def nbytes(self):
        return self.batch_buffers.nbytes

    def evaluate(self, points, wavenumbers):
        """Return S[k,point,panel], with resident geometry/tables reused."""
        self.set_points(points)
        ks = np.asarray(wavenumbers, dtype=float)
        if ks.ndim != 1 or not len(ks) or not np.isfinite(ks).all() or np.any(ks <= 0):
            raise ValueError('Prototype requires finite positive wavenumbers')
        output = np.empty((len(ks), self.np, self.nf), dtype=np.complex128)
        buffer = self.cl.Buffer(self.context, self.cl.mem_flags.WRITE_ONLY, output.nbytes)
        event = self._evaluate_into(ks, buffer)
        self.cl.enqueue_copy(self.queue, output, buffer, wait_for=[event]).wait()
        self.last_kernel_seconds = (event.profile.end-event.profile.start)*1e-9
        return output

    def _evaluate_into(self, ks, buffer):
        wavenumbers = self.batch_buffers.get('wavenumbers', ks.nbytes)
        self.cl.enqueue_copy(self.queue, wavenumbers, np.ascontiguousarray(ks)).wait()
        return self.potential_kernel(self.queue, (len(ks)*self.np*self.nf,), None,
            self.batch_buffers['points'], np.int32(self.np), np.int32(self.nf), self.centers, self.areas,
            np.int32(self.nq), self.quadrature, self.weights, wavenumbers, np.int32(len(ks)),
            np.int32(self.nr), np.int32(self.nz), self.rr, self.zz, self.table,
            self.pi, self.gamma, self.batch_buffers['rankine'], buffer)

    def apply_sources(self, points, wavenumbers, sources, *, _return_device=False):
        """Return only field responses; matrices stay on GPU in bounded batches."""
        points = np.asarray(points, dtype=float)
        if (points.ndim != 2 or points.shape[1] != 3 or not len(points)
                or not np.isfinite(points).all() or np.any(points[:, 2] != 0)):
            raise ValueError('GPU fields require finite nonempty (point,3) z=0 targets')
        ks = np.asarray(wavenumbers, dtype=float)
        source = np.ascontiguousarray(sources, dtype=np.complex128)
        if (ks.ndim != 1 or not len(ks) or not np.isfinite(ks).all() or np.any(ks <= 0)
                or source.ndim != 3 or source.shape[:2] != (len(ks), self.nf)
                or source.shape[2] < 1):
            raise ValueError('Sources must have (frequency,panel,mode) shape and positive wavenumbers')
        if source.flags.writeable or source is not self._source_values:
            if not np.isfinite(source).all():
                raise ValueError('Sources must be finite')
            self._source_values = None
            source_buffer = self.batch_buffers.get('sources', source.nbytes)
            self.cl.enqueue_copy(self.queue, source_buffer, source).wait()
            self._source_values = source
        budget = self.matrix_budget()
        if self.nf*16 > budget:
            raise RuntimeError('Source mesh exceeds the GPU matrix allocation budget')
        shape=(len(ks),len(points),source.shape[2])
        if _return_device:
            output=DeviceTransfer(self.batch_buffers.get('transfer', int(np.prod(shape))*16),shape,self.context)
        else:
            output = np.empty(shape, np.complex128)
        point_batch = max(1, budget//(self.nf*16))
        for p0 in range(0, len(points), point_batch):
            p1 = min(p0+point_batch, len(points))
            self.set_points(points[p0:p1])
            frequency_batch = max(1, budget//(self.np*self.nf*16))
            for f0 in range(0, len(ks), frequency_batch):
                f1 = min(f0+frequency_batch, len(ks))
                matrix = self.batch_buffers.get('matrix', (f1-f0)*self.np*self.nf*16)
                result = (output.buffer if _return_device else
                          self.batch_buffers.get('applied', (f1-f0)*self.np*source.shape[2]*16))
                event = self._evaluate_into(ks[f0:f1], matrix)
                applied = self.apply_kernel(self.queue, ((f1-f0)*self.np*source.shape[2],), None,
                    matrix, self.batch_buffers['sources'], np.int32(self.nf), np.int32(self.np),
                    np.int32(source.shape[2]), np.int32(f1-f0), np.int32(f0),
                    np.int32(len(points) if _return_device else self.np),
                    np.int32(p0 if _return_device else 0),np.int32(f0 if _return_device else 0),result,
                    wait_for=[event])
                if _return_device:
                    continue
                target = output[f0:f1, p0:p1]
                contiguous = target if target.flags.c_contiguous else np.empty(target.shape, np.complex128)
                self.cl.enqueue_copy(self.queue, contiguous, result, wait_for=[applied]).wait()
                if contiguous is not target:
                    target[:] = contiguous
        return output

    def apply_sources_device(self,points,wavenumbers,sources):
        return self.apply_sources(points,wavenumbers,sources,_return_device=True)

    apply_sources_transfer = apply_sources_device

    def incident_spectrum(self,train,xy):
        """Physical incident plane-wave sum; no change to phases or wave support."""
        xy=np.asarray(xy,dtype=float)
        if xy.ndim!=2 or xy.shape[1]!=2 or not np.isfinite(xy).all():
            raise ValueError('Incident targets must be finite (point,2) coordinates')
        coefficients=np.asarray(train.C)
        if coefficients.shape!=(len(train.theta),2,len(train.k)) or not np.isfinite(coefficients).all():
            raise ValueError('Incident coefficients must have (carrier,branch,frequency) shape')
        directions=np.c_[np.cos(train.theta),np.sin(train.theta)]
        packed=np.ascontiguousarray(coefficients.transpose(2,0,1),dtype=np.complex128)
        waves=self.cl.Buffer(self.context,self.cl.mem_flags.READ_ONLY|self.cl.mem_flags.COPY_HOST_PTR,
                             hostbuf=packed)
        output=np.empty((len(train.k),len(xy)),np.complex128)
        buffer=self.cl.Buffer(self.context,self.cl.mem_flags.WRITE_ONLY,output.nbytes)
        event=self.incident_kernel(self.queue,(output.size,),None,self.upload(xy),np.int32(len(xy)),
            self.upload(train.k),np.int32(len(train.k)),self.upload(directions),np.int32(len(directions)),
            waves,buffer)
        self.cl.enqueue_copy(self.queue,output,buffer,wait_for=[event]).wait()
        return output
