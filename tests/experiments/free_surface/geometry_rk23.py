"""RK23 trial rejection for the P4 sampled-cell Jacobian guard only.

Installed SciPy rk_step overwrites K[0] and all needed subsequent stage rows
on retry. Its _step_impl commits t/y/f only after a successful step; a failed
RHS therefore leaves the accepted state intact. All ordinary RK23 arithmetic
is delegated unchanged. Physical contact/rim failures and other errors escape.
"""
import numpy as np
from scipy.integrate import RK23

from .p4_geometry import InvalidGeometryError


class GeometryRK23(RK23):
    invalid_geometry_rejections = 0

    def _step_impl(self):
        original_fun = self.fun
        stage_time = None

        def trial_fun(time, state):
            nonlocal stage_time
            stage_time = time
            return original_fun(time, state)

        self.fun = trial_fun
        try:
            while True:
                minimum = 10*abs(np.nextafter(self.t, self.direction*np.inf)-self.t)
                proposal = min(max(self.h_abs, minimum), self.max_step, abs(self.t_bound-self.t))
                try:
                    return super()._step_impl()
                except InvalidGeometryError:
                    self.invalid_geometry_rejections += 1
                    # SciPy can reject for accuracy before the geometry failure
                    # without updating self.h_abs. The failing stage interval is
                    # <= its actual trial h, so half that interval is a safe
                    # conservative backoff, never a larger retry of that trial.
                    attempted = min(proposal, abs(stage_time-self.t))
                    reduced = attempted/2
                    if reduced < minimum or self.t+self.direction*reduced == self.t:
                        return False, ('Invalid sampled P4 trial geometry persists at the '
                                       'minimum representable RK23 step.')
                    self.h_abs = reduced
        finally:
            self.fun = original_fun
