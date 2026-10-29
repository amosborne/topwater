"""Build and check both distributions in isolated pip environments.

Run after building into ``~/.cache/topwater/dist``; use ``--directory`` to override.
Temporary environments and the rebuilt wheel are removed when the check finishes.
"""

import argparse
from pathlib import Path
import subprocess
import tempfile
import tomllib
import venv


SMOKE_CHECK = """
from importlib.metadata import metadata
from importlib import resources, util
from pathlib import Path
import sys
import numpy
import scipy
import topwater
from topwater.sea import states as sea_states
from topwater.boat.generation import generate_boat
from topwater.cli import parser
from topwater.workflow_data import Boat, Sea, Result

info = metadata('topwater')
assert info['Name'] == 'topwater'
assert info['Version'] == sys.argv[1]
requirements = info.get_all('Requires-Dist', [])
assert any(item.startswith('numpy>=') for item in requirements)
assert any(item.startswith('scipy>=') for item in requirements)
assert Path(topwater.__file__).resolve().is_relative_to(Path.cwd() / 'env')
assert resources.files('topwater.sea').joinpath('data/socal_sea_states.json').is_file()
assert resources.files('topwater.solve.wecsim').joinpath('data/run_workflow_case.m').is_file()
assert resources.files('topwater.sea').joinpath('data/stronger-sea-observation.npz').is_file()
for name in ('backend.py', 'frames.py', 'workflow.py', 'geometry.py', 'scene.py', 'mesh.py', 'hud.py', 'coolwarm.npy', 'NOTICE.txt', 'PALETTE-LICENSE.txt'):
    assert resources.files('topwater.render').joinpath(name).is_file()
import skia
from PIL import Image
assert Image.new('RGBA', (2, 2)).size == (2, 2)
assert util.find_spec('wgpu') is None
from topwater.render.backend import DrawingBackend
drawing = DrawingBackend((32, 32), [-1., -1.], [1., 1.])
try:
    vertices = numpy.array([[-1, -1, 1, 0, 0, 1], [1, -1, 1, 0, 0, 1], [0, 1, 1, 0, 0, 1]], numpy.float32)
    pixels = drawing.draw(vertices, numpy.array([[0, 1, 2]], numpy.uint32), [], [])
    assert pixels.shape == (32, 32, 4)
    assert pixels[16, 16, 0] > 200 and pixels[16, 16, 1] < 20
finally:
    drawing.close()
for owner in ('cpu', 'gpu'):
    assert resources.files('topwater.render.' + owner).joinpath('backend.py').is_file()
assert sea_states.catalog()
defaults = parser().parse_args(['solve'])
assert defaults.solver == 'native' and defaults.dwell == 100. and defaults.time_step_s == .01
from topwater.solve import SolveOptions
from topwater.solve.settings import VIDEO_SIZE_PX, SURFACE_PIXELS_PER_SAMPLE
from topwater.solve import LoadEvaluation
assert util.find_spec('topwater.native') is None
assert SolveOptions().render_every == 10
assert VIDEO_SIZE_PX == (800, 800) and SURFACE_PIXELS_PER_SAMPLE == 10
assert SolveOptions().surface_interval_s == .1 and SolveOptions().packet_tolerance_m == .01
assert util.find_spec('pyopencl') is None
assert util.find_spec('matlab') is None
from topwater.solve.cpu.backend import CPUFieldBackend
assert 'pyopencl' not in sys.modules
for name in ('fields.py', 'source_cache.py', 'transfer.py', 'settings.py', 'sampling.py', 'surface.py', 'hydrodynamics.py', 'diffraction.py', 'validation.py', 'progress.py', 'types.py', 'fourier.py', 'recovery.py', 'analysis.py', 'evidence.py'):
    assert resources.files('topwater.solve').joinpath(name).is_file()
for name in ('potential.cl', 'fourier.cl', 'LICENSE'):
    assert resources.files('topwater.solve.gpu').joinpath(name).is_file()
for name in ('backend.py', 'potential.py', 'fourier.py'):
    assert resources.files('topwater.solve.cpu').joinpath(name).is_file()
boat = generate_boat(7.3152, 2.7432, 20, 2500)
assert boat.mass_kg == 2500
for retired in ('bem', 'boundary', 'inertia', 'rankine', 'fluid_boundary', 'surface', 'native', 'result_render', 'visualization', '_visualization_2d', '_visualization_geometry'):
    assert util.find_spec('topwater.' + retired) is None
print(f"Installed {info['Name']} {info['Version']} from distribution")
"""


def only_file(directory, pattern):
    matches = list(directory.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one {pattern} in {directory}, found {len(matches)}")
    return matches[0]


def check_wheel(wheel, directory, version):
    directory.mkdir()
    environment = directory / "env"
    venv.EnvBuilder(with_pip=True).create(environment)
    python = environment / (
        "Scripts/python.exe" if (environment / "Scripts").exists() else "bin/python"
    )
    subprocess.run(
        [str(python), "-I", "-m", "pip", "install", str(wheel)],
        cwd=directory,
        check=True,
    )
    subprocess.run(
        [str(python), "-I", "-c", SMOKE_CHECK, version], cwd=directory, check=True
    )


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path,
                        default=Path.home() / '.cache/topwater/dist',
                        help='Directory containing the built wheel and source archive')
    args = parser.parse_args()
    version = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))[
        "project"
    ]["version"]
    wheel = only_file(args.directory, "*.whl")
    source = only_file(args.directory, "*.tar.gz")
    with tempfile.TemporaryDirectory(prefix="topwater-packaging-") as temporary:
        directory = Path(temporary).resolve()
        check_wheel(wheel, directory / "wheel-check", version)
        rebuilt = directory / "rebuilt"
        subprocess.run(
            ["uv", "build", str(source), "--wheel", "--out-dir", str(rebuilt)],
            cwd=directory,
            check=True,
        )
        check_wheel(only_file(rebuilt, "*.whl"), directory / "source-check", version)


if __name__ == "__main__":
    main()
