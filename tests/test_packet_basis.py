"""Identical packet numerics, bounded memory, and per-grid invalidation."""
import numpy as np

from topwater.solve.packet_basis import PacketBasis
from topwater.solve.packets import minimum_cutoff_packet


def test_reused_basis_preserves_cutoff_and_every_diagnostic():
    cache=PacketBasis()
    for frequency,phase,amplitude in ((.25,.4,.2),(.43,-.8,1.2),(.25,1.,.8)):
        plain=minimum_cutoff_packet(frequency,phase,1.,1.,4.,amplitude_m=amplitude)
        reused=minimum_cutoff_packet(frequency,phase,1.,1.,4.,amplitude_m=amplitude,_harmonics=cache)
        assert plain.cutoff_hz==reused.cutoff_hz
        np.testing.assert_array_equal(plain.coefficients,reused.coefficients)
        np.testing.assert_array_equal(plain.elevation_m,reused.elevation_m)
        for name in plain.diagnostics:
            np.testing.assert_array_equal(plain.diagnostics[name],reused.diagnostics[name])


def test_budget_and_grid_invalidation():
    cache=PacketBasis(maximum_bytes=160)
    time=np.arange(10)*.1
    row=cache.prepare(time,1.)(1)
    assert cache.prepare(time.copy(),1.)(1) is row
    cache.harmonic(2)
    assert sum(value.nbytes for value in cache.rows.values())<=160
    shifted=time+.01
    changed=cache.prepare(shifted,1.)(1)
    np.testing.assert_array_equal(changed,np.exp(2j*np.pi*shifted))
    assert changed is not row
