# Dependency licenses

Topwater's own source is MIT licensed. Dependencies retain their own licenses;
consult the notices distributed with the exact versions in `uv.lock`.

| Direct dependency | Role | Upstream license |
| --- | --- | --- |
| [NumPy](https://github.com/numpy/numpy/blob/main/LICENSE.txt) | Numerical arrays; runtime | BSD-3-Clause |
| [SciPy](https://github.com/scipy/scipy/blob/main/LICENSE.txt) | Numerical routines; runtime | BSD-3-Clause |
| [pytest](https://github.com/pytest-dev/pytest/blob/main/LICENSE) | Development tests | MIT |
| [Hatchling](https://github.com/pypa/hatch/blob/master/LICENSE.txt) | PEP 517 build backend | MIT |

[uv](https://github.com/astral-sh/uv#license) is the development/build front end,
licensed under MIT or Apache-2.0. It is not a runtime dependency. Binary NumPy and
SciPy wheels can include additional libraries and accompanying license notices.
No third-party implementation source is copied into this repository.
