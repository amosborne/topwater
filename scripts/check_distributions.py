"""Build and check both distributions in isolated pip environments.

Run after ``uv build`` with ``uv run --locked python scripts/check_distributions.py``.
Temporary environments and the rebuilt wheel are removed when the check finishes.
"""

from pathlib import Path
import subprocess
import tempfile
import tomllib
import venv


SMOKE_CHECK = """
from importlib.metadata import metadata
from pathlib import Path
import sys
import numpy
import scipy
import topwater
from topwater.bem import assemble_dirichlet_halfspace
from topwater.hydrostatics import equilibrate_heave
from topwater.inertia import instantaneous_added_mass

info = metadata('topwater')
assert info['Name'] == 'topwater'
assert info['Version'] == sys.argv[1]
requirements = info.get_all('Requires-Dist', [])
assert any(item.startswith('numpy>=') for item in requirements)
assert any(item.startswith('scipy>=') for item in requirements)
assert Path(topwater.__file__).resolve().is_relative_to(Path.cwd() / 'env')
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
    version = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))[
        "project"
    ]["version"]
    wheel = only_file(root / "dist", "*.whl")
    source = only_file(root / "dist", "*.tar.gz")
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
