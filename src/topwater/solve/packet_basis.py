"""Per-case, bounded reuse of identical Fourier basis vectors.

Only numerical exponentials are cached; packet targets, cutoffs, bounds and
diagnostics are still computed independently for every carrier.
"""
import numpy as np


class PacketBasis:
    def __init__(self, maximum_bytes=64*1024**2):
        self.maximum_bytes = maximum_bytes
        self.time_s = None
        self.period_s = None
        self.rows = {}

    def prepare(self, time_s, period_s):
        if (self.time_s is None or self.period_s != period_s
                or not np.array_equal(time_s, self.time_s)):
            self.time_s = np.array(time_s, copy=True)
            self.period_s = period_s
            self.rows = {}
        return self.harmonic

    def harmonic(self, index):
        if index in self.rows:
            return self.rows[index]
        # Identical expression and order to the uncached packet algorithm.
        row = np.exp(2j*np.pi*index*self.time_s/self.period_s)
        if (len(self.rows)+1)*row.nbytes <= self.maximum_bytes:
            row.flags.writeable = False
            self.rows[index] = row
        return row
