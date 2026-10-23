# Immutable integer plane ownership

The exact integer cold path previously reconstructed and compared every plane
tuple on each owner access. A root whole-case profile found762228 checks taking
12.61s cumulatively under instrumentation. That profile motivated this change;
its time is not an unprofiled speed prediction.

On first integer use, PlanePoints validates its planes and adopts the integer
owner's immutable tuple of immutable rows. Later validation uses identity.
Replacing the registry planes still performs a coefficient comparison: equal
mutable replacements are re-frozen, changed replacements reject even on cached
intersection/sign queries before cache/counter changes. Adding a plane after
witness-only integer initialization also rejects. Fraction mode is unchanged.

Twelve targeted tests passed in5.92s on Python3.12.8. An initial test invocation
had sandbox temporary-directory setup errors; the successful run used a fresh
workspace basetemp. No physical test failed.

The single bounded CPU audit used old/new/new/old at generic t0 and the saved
actual.964s state. The old registry class was loaded from git73bbdad and retained
its original list of planes. Both modes used the same current integer helper;
the old mode never satisfies its identity fastpath and therefore executes the
original deep comparison. Both original helper and registry sources were
archived. This is a no-freeze integer control, not the slower Fraction oracle.

| Input | Old seconds | New seconds | Old/new median |
|---|---|---|---|
| Generic t0 | .75858,.81889 | .67705,.69536 |1.1494|
| Actual.964s | .83997,.88535 | .78951,.74295 |1.1258|

All eight calls matched exact saved arrays, templates, ordered polygons,
intersection registration, incidence/relevance and sign caches. Registry plane
containers were normalized from list/tuple for comparison; coefficients were
unchanged. Arrays also matched separately SHA-owned first-stage frozen outputs.
No fluid solve, force evaluation or moving trajectory was run. Observed cold
gains do not establish a universal factor or full-case throughput benefit.

The raw report is [boat-integer-freeze-cold.json](boat-integer-freeze-cold.json).
It records inputs, versions, exact producer hashes, caps and storage counts.
Ignored source archives and per-call outputs are retained under the isolated
worktree's artifacts/integer-freeze-cold-audit; a copy of the report and seven
current plus two original source owners is in root artifacts. No numerical
producer was edited after this audit.
