"""Immediate allocations and reusable, named per-solve buffers."""
import sys


class GPUAllocationError(RuntimeError):
    """A device-buffer allocation failed before field computation."""


class ReusableBuffers:
    def __init__(self, queue):
        import pyopencl as cl
        from pyopencl.tools import ImmediateAllocator
        self.cl, self.queue = cl, queue
        self.allocator = ImmediateAllocator(queue)
        self.values = {}

    def __getitem__(self, name):
        return self.values[name]

    def get(self, name, size):
        current = self.values.get(name)
        if current is not None and current.size >= size:
            return current
        if current is not None:
            self.queue.finish()
            self.values.pop(name).release()
        if size > self.queue.device.max_mem_alloc_size:
            raise GPUAllocationError(f'{name}: {size} bytes exceeds the device allocation limit')
        try:
            value = self.allocator(size)
        except self.cl.Error as error:
            if error.code not in (self.cl.status_code.MEM_OBJECT_ALLOCATION_FAILURE,
                                  self.cl.status_code.OUT_OF_RESOURCES,
                                  self.cl.status_code.OUT_OF_HOST_MEMORY):
                raise
            raise GPUAllocationError(f'{name}: could not allocate {size} bytes: {error}') from error
        self.values[name] = value
        return value

    @property
    def nbytes(self):
        return sum(value.size for value in self.values.values())

    def clear(self):
        try:
            self.queue.finish()
        finally:
            primary = sys.exception()
            values, self.values = self.values, {}
            errors = []
            for value in values.values():
                try:
                    value.release()
                except Exception as error:
                    errors.append(error)
            if errors:
                if primary is not None:
                    for error in errors:
                        primary.add_note(f'Buffer release also failed: {error}')
                else:
                    for error in errors[1:]:
                        errors[0].add_note(f'Buffer release also failed: {error}')
                    raise errors[0]
