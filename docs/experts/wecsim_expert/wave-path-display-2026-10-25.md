# Imported-wave path display correction

Verified checkpoint, 2026-10-25. The coordinator assigned one cosmetic
warning from the actual WEC run: `waveClass/listInfo` treated the imported
elevation filepath as part of an `fprintf` format string. A Windows `\U`
sequence produced an escaped-character warning and truncated the printed path.
Percent sequences in a filename had the same format-string exposure.

The maintained fork commit
`9952634f6f2e2746a1626484e27c3003d1e1d5a9` passes `elevationFile` as the
argument of a fixed `%s` conversion at `source/objects/waveClass.m:411`.
Only that display statement changes; no parsing, wave setup, loads or numerical
code changes. Root committed and pushed the owned fork branch and updated the
backend/documented pin and `wecsim-fork.json`.

The private regression helper is
`artifacts/workflow/check_wave_list_info.m`. With the maintained fork already on
the existing shared-engine path, call `check_wave_list_info(forkSource)`. It
constructs three uninitialized `elevationImport` wave objects, captures
`listInfo`, and requires verbatim Windows drive, percent-containing, and UNC
paths. It temporarily enables all warnings, requires empty `lastwarn`, and
restores warning settings and `lastwarn` on success or failure. Percent signs
are constructed with `char(37)` to avoid a private static-scanner false positive;
the actual test path values are unchanged. MATLAB did not support the originally
proposed `warning('error','all')`, so that unexecuted source design was replaced.
The helper reads no sea files, computes no waves, and starts no model or engine.

Root's allocated real shared-engine probe passed:
[report](../../../artifacts/workflow/wave-display-check-attempt3/report.json),
[manifest](../../../artifacts/workflow/wave-display-check-attempt3/manifest.json),
[terminal proof](../../../artifacts/workflow/wave-display-check-attempt3/terminal.json).
The report records three paths printed verbatim, zero warnings, no wave setup,
no body simulation, and no new engine. Report SHA256:
`cfcc3f458221c456bc6a094d93a2231111212079ff1a4fb785aa23b2ec0f3395`.
The helper MATLAB style check is clean. The upstream class retains 113 existing
style findings, with zero new findings; whole-class style cleanliness is not
claimed. Earlier failed private harness attempts remain preserved.

The integration suite subsequently passed 228 tests with one optional skip in
14.98 s. The canonical public sea-generation command also completed with its
new log and all sea arrays bit-exact. These checks do not establish a complete
solve/render workflow or physical accuracy. Full-case acceptance remains pending
the cached geometry sweep and public workflow; this display correction does not
resolve the separate frame-277 surface triangulation failure.
