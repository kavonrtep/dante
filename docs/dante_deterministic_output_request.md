# Feature request (DANTE): deterministic domain output

**Status:** requested — fix belongs in `DANTE`, not CARP. **Priority: low.**
**Affects:** `DANTE.gff3` (protein_domain records).
**CARP tracking:** part of `fix/deterministic-library-clustering`.

## Summary

`DANTE.gff3` is **not byte-reproducible** run-to-run at a fixed version
(DANTE 0.2.11, same genome, same host). The annotation is *almost* stable, but a
small fraction of records differ, which is enough to seed downstream
order-sensitive clustering (see the companion dante_ltr request) into large
library churn.

## Evidence (measured in CARP)

Two runs of DANTE 0.2.11 on the same input:

- **Feature order identical**, and across 2.74 M filtered domains **0
  `Final_Classification` differences** — so the primary call is deterministic.
- A `cmp` still differs: the per-domain **`Region_Hits_Classifications`** list is
  reordered and its top hit occasionally changes between near-tie REXdb lineages
  (e.g. `RH|…|Ty1/copia|Lyco[195bp]` vs `RH|…|Ty1/copia|SIRE[207bp]` at the same
  locus).
- Repeated single-sequence reruns (D1…D5): identical feature count (28065),
  identical order, **0 `Final_Classification` changes**, but ~0.1 % of features
  (≈27/28065) have slightly jittered boundaries.

DANTE runs multithreaded; the variability looks like thread-order-dependent
resolution of near-equal-scoring hits.

## Requested fix (in DANTE)

Make the output byte-deterministic at fixed input+version:

1. Break ties in best-hit / `Region_Hits_Classifications` selection with a fixed
   rule (e.g. score desc, then alignment length desc, then lineage name), and
   emit `Region_Hits_Classifications` in a sorted, canonical order.
2. Ensure domain-boundary determination does not depend on thread scheduling.

## Impact / priority

Low on its own: `Final_Classification` is already deterministic, so DANTE's
*classification* output is stable. The residual boundary jitter only matters
because CARP's clustering was order-sensitive — which CARP now fixes on its own
side by sorting clustering inputs. This request removes the seed at the source
and gives byte-reproducible DANTE output, but is not required for CARP library
reproducibility once the sorting (CARP-side) and the dante_ltr fix are in place.

## Resolution (DANTE 0.2.12)

Requested fix **#1 — done.** The run-to-run non-determinism was diagnosed as
`list(set(...))` iteration order in `group_annot_regs` (`dante.py`), which
varies with `PYTHONHASHSEED` between processes and reshuffled
`Region_Hits_Classifications`. Fixed by:

- `group_annot_regs` now emits `unique_annotations` and `domain_type` in
  canonical (`sorted`) order, with `ann_pos_counts` kept paired.
- `best_score` gained an explicit deterministic tie-break: score desc,
  alignment length desc, lineage name, then DB sequence name.

Verified: two runs of the same input with different `PYTHONHASHSEED` are now
byte-identical (baseline `HEAD` differed); guarded by `tests.sh` +
`.github/workflows/tests.yml`.

Requested fix **#2 — not reproduced / not done.** Static analysis showed DANTE
runs `lastal` single-threaded (`-p` is the scoring matrix, not `-P` threads);
parallelism is a process `Pool` over genome chunks whose output is globally
sorted + de-duplicated, so the final bytes do not depend on thread scheduling.
The residual ~0.1 % boundary jitter is therefore most likely `lastal`
alignment tie-breaking or overlap-chunk duplicate resolution, not Python thread
order. It needs an empirical repro before any change and was left out of this
(low-priority) fix.
