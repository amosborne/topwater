"""Versioned, exact numeric radiation-source reuse between solves."""
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import sys
import tempfile
from zipfile import BadZipFile

import numpy as np
from .types import MODES


class SourceCache:
    def __init__(self, directory=None):
        self.directory = Path(directory or os.environ.get('TOPWATER_SOURCE_CACHE_DIR')
                              or Path.home()/'.cache/topwater/radiation-sources-v1')

    def path(self, hydro, omega):
        mesh = hydro.body.mesh_including_lid.merged()
        green = hydro.solver.engine.green_function
        digest = hashlib.sha256()
        metadata = dict(schema=1, capytaine=version('capytaine'), numpy=np.__version__,
                        scipy=version('scipy'), density=hydro.density, gravity=hydro.gravity,
                        depth=hydro.depth, engine=repr(hydro.solver.engine), green=repr(green))
        digest.update(json.dumps(metadata, sort_keys=True).encode())
        arrays = [omega, mesh.vertices, mesh.faces, mesh.faces_centers, mesh.faces_normals,
                  mesh.faces_areas, mesh.faces_radiuses, *mesh.quadrature_points,
                  green.tabulated_r_range, green.tabulated_z_range, green.tabulated_integrals]
        for name in MODES:
            arrays.append(hydro.body.dofs[name].evaluate_motion(hydro.body.mesh))
        for array in arrays:
            values = np.ascontiguousarray(array)
            if values.dtype.hasobject:
                raise ValueError('Source cache identity requires numerical arrays')
            digest.update(str((values.shape,values.dtype.str)).encode())
            digest.update(values.tobytes())
        return self.directory/(digest.hexdigest()+'.npz')

    def load(self, path, shape):
        try:
            with np.load(path,allow_pickle=False) as data:
                values=data['sources']
                if values.shape!=shape or values.dtype!=np.complex128 or not np.isfinite(values).all():
                    raise ValueError('invalid source-bank arrays')
            values.flags.writeable=False
            print(f'SOURCE_CACHE hit key={path.stem}',file=sys.stderr)
            return values
        except FileNotFoundError:
            return None
        except (OSError,ValueError,KeyError,BadZipFile,EOFError) as error:
            print(f'SOURCE_CACHE rebuild key={path.stem} reason={error}',file=sys.stderr)
            return None

    def save(self,path,sources):
        staged=None
        try:
            path.parent.mkdir(parents=True,exist_ok=True)
            fd,name=tempfile.mkstemp(prefix='.sources-',suffix='.npz',dir=path.parent)
            staged=Path(name)
            with os.fdopen(fd,'wb') as stream:
                np.savez(stream,sources=sources)
            os.replace(staged,path)
            return True
        except OSError as error:
            print(f'SOURCE_CACHE write skipped: {error}',file=sys.stderr)
            return False
        finally:
            if staged is not None:
                try:
                    staged.unlink(missing_ok=True)
                except OSError as error:
                    print(f'SOURCE_CACHE staging cleanup skipped at {staged}: {error}',file=sys.stderr)
