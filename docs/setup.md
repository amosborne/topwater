# Setup

Python 3.12 or newer and [uv](https://docs.astral.sh/uv/) are required.
Clone the repository and fetch its Git LFS files to use the published example:

```powershell
git lfs install
git lfs pull
uv sync --locked
uv run topwater --help
```

The base installation supports input/output, sea generation, plotting and video
rendering. The renderer uses the FFmpeg executable supplied by imageio-ffmpeg.

## Simulation dependencies

The supported simulation setup is MATLAB R2025b, matlabengine 25.2.2,
Capytaine 3.0.0 and licensed Simulink, Simscape and Simscape Multibody.
Git must be available to verify the WEC-Sim revision.

```powershell
uv sync --locked --group wecsim
git clone https://github.com/amosborne/WEC-Sim.git "$HOME/.cache/topwater/wec-sim"
git -C "$HOME/.cache/topwater/wec-sim" checkout 6b774228e1f6b31dee473e626204b0ef65e5ec38
```

The maintained fork is based on WEC-Sim v7.1.0. Its exact revision is checked
before use; modified source or another revision is rejected. Topwater's default
checkout location is `~/.cache/topwater/wec-sim`. To use a different checkout:

```powershell
$env:TOPWATER_WECSIM_PATH = 'C:/path/to/WEC-Sim'
```

Topwater attaches to a shared MATLAB engine named `CAPYTAINE_REFERENCE`.
On Windows, start it with:

```powershell
./scripts/start_matlab.ps1
```

The helper launches MATLAB hidden, sets separate temporary/cache/code-generation
directories under `~/.cache/topwater/matlab/`, and shares the named engine.
It refuses to launch while MATLAB is already running so an existing session can
be reconciled first. Alternatively, configure a MATLAB session yourself and run:

```matlab
matlab.engine.shareEngine('CAPYTAINE_REFERENCE')
```

Keep MATLAB's temporary, Simulink cache and code-generation directories as
separate siblings. In particular, TMP/TEMP must not live inside the initial
Simulink CacheFolder. The adapter restores the caller's engine settings and
reports cleanup failures instead of publishing a result from an incomplete run.
It attaches to the engine; the caller owns starting and stopping MATLAB.
Do not overlap solves on the same shared engine.

## Files generated during use

Published files go into the case directory selected on the command line.
Intermediate BEMIO, MATLAB, Simulink and surface products use OS temporary storage.
Failed runs can retain a temporary recovery location and print its path.
Reusable dependencies/caches live outside the checkout under `~/.cache/topwater/`.
The solver reuses Capytaine's numerical Green-function tabulations under
`~/.cache/topwater/capytaine-3.0.0`; boat coefficients and motion are recomputed.
There is no repository-local artifacts directory.

## Local development

```powershell
uv sync --locked --group dev
uv run --locked pytest -p no:cacheprovider
uv build --out-dir "$HOME/.cache/topwater/dist"
uv run --locked python scripts/check_distributions.py --directory "$HOME/.cache/topwater/dist"
```

Tests exercise implementation behavior; they are not a physical qualification
certificate. MATLAB contract tests under `tests/matlab` are optional native
regressions and are not invoked by ordinary pytest. No CI workflow is configured.

See [performance](performance.md) for the simple-boat profiling harness,
measured speedups and remaining opportunities.
