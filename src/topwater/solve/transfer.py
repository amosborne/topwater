"""Borrowed per-batch transfer storage with shared lazy scaling and DC row."""
from dataclasses import dataclass, replace


@dataclass(frozen=True)
class Transfer:
    buffer: object
    shape: tuple
    context: object = None
    scale: float = 1.
    leading_zero: bool = False

    def scaled(self, value):
        return replace(self, scale=self.scale*value)

    def with_zero(self):
        return replace(self, shape=(self.shape[0]+1, *self.shape[1:]), leading_zero=True)


def prepare_host_transfer(transfer, nw):
    from topwater.solve.validation import _finite
    if isinstance(transfer, Transfer):
        if transfer.context is not None or transfer.shape[0] != nw or transfer.shape[-1] != 6:
            raise ValueError('CPU transfer must use host storage and the requested frequency/mode layout')
        _finite(transfer.buffer, (nw-int(transfer.leading_zero), transfer.shape[1], 6),
                'velocity_transfer', complex_value=True)
        return transfer
    return _finite(transfer, (nw, None, 6), 'velocity_transfer', complex_value=True)
