"""One owner for reconstruction caching, batching and source/result coupling."""
from functools import wraps
import sys
import numpy as np
from topwater.solve.types import MODES
from topwater.solve.progress import _summarize_capytaine_warnings
from topwater.solve.source_cache import SourceCache
from topwater.solve.fourier import fft_dimensions


def _field_errors(function):
    @wraps(function)
    def call(*args, **kwargs):
        try:
            return function(*args, **kwargs)
        except Exception as error:
            cl = getattr(args[0].evaluator, 'cl', None)
            if cl is not None and isinstance(error, cl.Error):
                raise RuntimeError(f'{args[0].device.upper()} {function.__name__} failed: {error}') from error
            raise
    return call


class FieldBackend:
    def __init__(self, hydro, evaluator, fourier, *, cache_sources=True):
        self.hydro, self.evaluator, self.fourier = hydro, evaluator, fourier
        self._sources = self._frequencies = None
        self.source_cache = SourceCache() if cache_sources else None

    @_field_errors
    def prepare_radiation(self, stamp, output_time, omega, settings, point_count):
        """Reserve the case's real workspaces before computing spatial chunks."""
        n, _ = fft_dimensions(len(stamp), stamp[1]-stamp[0], settings)
        points = min(settings.chunk_points, point_count)
        while True:
            try:
                self.fourier.reserve(n, points, len(output_time))
                self.evaluator.reserve_radiation(points, len(omega)-1)
                print(f'{self.device.upper()}_WORKSPACE fft_samples={n} fft_length={self.fourier.m} '
                      f'batch_points={points} reusable_bytes='
                      f'{self.fourier.nbytes+self.evaluator.nbytes}', file=sys.stderr)
                return points
            except self.evaluator.allocation_error as error:
                self.fourier.reset_buffers()
                self.evaluator.reset_batch_buffers()
                if points == 1:
                    raise RuntimeError(f'{self.device.upper()} workspace cannot fit a one-point batch: {error}') from error
                reduced = max(1, points//2)
                print(f'{self.device.upper()}_WORKSPACE allocation failed at batch_points={points}; '
                      f'retrying batch_points={reduced}: {error}', file=sys.stderr)
                points = reduced

    @_field_errors
    def finish_radiation(self):
        """Release stage-owned workspaces before diffraction and incident fields."""
        try:
            self.fourier.reset_buffers()
        finally:
            primary = sys.exception()
            try:
                self.evaluator.reset_batch_buffers()
            except Exception as error:
                if primary is None:
                    raise
                primary.add_note(f'{self.device.upper()} potential workspace cleanup also failed: {error}')

    @_summarize_capytaine_warnings
    @_field_errors
    def radiation_transfer(self, hydro, xy, omega):
        import capytaine as cpt
        if self._frequencies is None or not np.array_equal(omega, self._frequencies):
            shape = (len(omega), hydro.body.mesh_including_lid.nb_faces, 6)
            path = None if self.source_cache is None else self.source_cache.path(hydro, omega)
            sources = None if path is None else self.source_cache.load(path, shape)
            if sources is None:
                sources = np.empty(shape, complex)
                for iw, w in enumerate(omega):
                    for j, mode in enumerate(MODES):
                        key = (float(w), mode)
                        if key not in hydro.radiation_results:
                            hydro.radiation_results[key] = hydro.solver.solve(cpt.RadiationProblem(
                                body=hydro.body, rho=hydro.density, g=hydro.gravity,
                                water_depth=hydro.depth, omega=float(w), radiating_dof=mode), keep_details=True)
                        sources[iw, :, j] = hydro.radiation_results[key].sources
                if path is not None:
                    if self.source_cache.save(path,sources):
                        print(f'SOURCE_CACHE stored key={path.stem}',file=sys.stderr)
            sources.flags.writeable = False
            self._sources = sources
            self._frequencies = np.array(omega, copy=True)
        response = self.evaluator.apply_sources_transfer(np.c_[xy, np.zeros(len(xy))],
                                                np.asarray(omega)**2/hydro.gravity, self._sources)
        return response.scaled(-1/hydro.gravity)

    def prepend_zero(self,transfer):
        return transfer.with_zero()

    @_field_errors
    def free_surface_elevation(self, xy, result):
        if result.forward_speed != 0:
            raise ValueError('Surface fields require zero forward speed')
        response = self.evaluator.apply_sources(np.c_[xy, np.zeros(len(xy))],
            np.array([result.encounter_wavenumber]), np.asarray(result.sources)[None, :, None])[0, :, 0]
        return -1/result.g*(-1j*result.omega)*response

    @_field_errors
    def incident_spectrum(self, train, xy):
        return self.evaluator.incident_spectrum(train, xy)

    @_field_errors
    def close(self):
        try:
            self.fourier.close()
        finally:
            primary = sys.exception()
            try:
                self.evaluator.close()
            except Exception as error:
                if primary is None:
                    raise
                primary.add_note(f'{self.device.upper()} evaluator cleanup also failed: {error}')
