"""Surface output samples selected from the accepted motion clock."""
import numpy as np
from topwater.solve.settings import DEFAULTS


def surface_sample_times(time_s, *, every=DEFAULTS.render_every):
    """Save regular native samples plus the exact, possibly partial endpoint."""
    time = np.asarray(time_s)
    if not isinstance(every, int) or isinstance(every, bool) or every < 1:
        raise ValueError('Surface sample stride must be a positive integer')
    indices = np.arange(0, len(time), every)
    if indices[-1] != len(time) - 1:
        indices = np.r_[indices, len(time) - 1]
    return time[indices]
