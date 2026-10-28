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
from topwater import sea_states
from topwater.cli import parser
from topwater.workflow_data import Boat, Sea, Result

info = metadata('topwater')
assert info['Name'] == 'topwater'
assert info['Version'] == sys.argv[1]
requirements = info.get_all('Requires-Dist', [])
assert any(item.startswith('numpy>=') for item in requirements)
assert any(item.startswith('scipy>=') for item in requirements)
assert Path(topwater.__file__).resolve().is_relative_to(Path.cwd() / 'env')
assert resources.files('topwater').joinpath('data/socal_sea_states.json').is_file()
assert resources.files('topwater').joinpath('data/wecsim/run_workflow_case.m').is_file()
assert resources.files('topwater').joinpath('data/stronger-sea-observation.npz').is_file()
assert sea_states.catalog()
assert parser().parse_args(['solve', '100']).solver == 'wecsim'
for retired in ('bem', 'boundary', 'inertia', 'rankine', 'fluid_boundary', 'surface'):
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
