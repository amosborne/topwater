# Capytaine-derived potential formulas

`potential.cl` and `../cpu/potential.py` translate the potential portions of Capytaine 3.0.0's
`libs/libDelhommeau/src/Green_Rankine.f90`, `Green_wave.f90`,
`Delhommeau_integrals.f90`, and `matrices.f90`.

Copyright 2026 Capytaine developers. Licensed under the Apache License, Version
2.0: <https://www.apache.org/licenses/LICENSE-2.0>.

The source distribution is available from <https://pypi.org/project/capytaine/3.0.0/>.
This is an independently implemented optional backend; the installed
Capytaine distribution is not modified. Source formulas have been translated
to double-precision OpenCL C and Numba CPU kernels, derivatives removed, frequency-independent terms
reused, and the supported domain restricted to the current surface workflow.
