"""Exact coherent forcing at the finite packet frequencies and directions."""
from dataclasses import dataclass

import numpy as np
from threadpoolctl import threadpool_limits

from topwater.solve.diffraction import packet_diffraction


@dataclass(frozen=True)
class Excitation:
    omega: np.ndarray
    coefficients: np.ndarray

    @classmethod
    def from_case(cls, setup, hydro):
        from topwater.solve.progress import _StageProgress
        train = setup.train
        response = packet_diffraction(hydro, train)
        coefficients = np.zeros((len(train.omega), 6), complex)
        progress = _StageProgress('native packet excitation', len(train.omega)-1, 'packet frequencies')
        with threadpool_limits(1):
            for i in range(1, len(train.omega)):
                coefficients[i] = response.force(i)
                progress.update(i)
        progress.finish()
        result = cls(np.asarray(train.omega), coefficients)
        scale = np.sum(abs(coefficients), axis=0)
        q0 = np.array([0., 0., -setup.world_vertical_offset_m, 0., 0., 0.])
        net = result.at(0.)-hydro.K@q0
        weight = setup.body.mass*hydro.gravity
        length = np.ptp(setup.body.mesh.vertices, axis=0).max()
        if (np.linalg.norm(net[:3]) > .001*weight
                or np.linalg.norm(net[3:]) > .001*weight*length
                or np.any(abs(net) > .001*scale)):
            raise ValueError('Native initial excitation plus restoring exceeds calm startup limits')
        return result

    def at(self, time_s):
        return (np.exp(-1j*self.omega*time_s)@self.coefficients).real
