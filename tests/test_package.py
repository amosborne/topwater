"""Check the installed development package against its declared metadata."""

from importlib.metadata import metadata
from pathlib import Path
import tomllib

import topwater


def test_package_metadata():
    project = tomllib.loads(
        (Path(__file__).resolve().parents[1] / "pyproject.toml").read_text(
            encoding="utf-8"
        )
    )["project"]
    installed = metadata("topwater")
    assert installed["Name"] == project["name"]
    assert installed["Version"] == project["version"]
    assert installed["Requires-Python"] == project["requires-python"]
    assert installed["License-Expression"] == project["license"]
    assert set(installed.get_all("Requires-Dist", [])) == set(project["dependencies"])
    assert topwater.__file__ is not None
