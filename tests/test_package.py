"""Check the installed development package against its declared metadata."""

from importlib.metadata import metadata
from pathlib import Path
import tomllib
from packaging.requirements import Requirement

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
    declared = installed.get_all("Requires-Dist", [])
    core = {item for item in declared if 'extra ==' not in item}
    assert {Requirement(item) for item in core} == {Requirement(item) for item in project["dependencies"]}
    assert set(installed.get_all("Provides-Extra", [])) == set(project.get("optional-dependencies", {}))
    assert topwater.__file__ is not None
