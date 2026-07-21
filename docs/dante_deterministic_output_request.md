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

Requested fix **#2 — disproved on production data; nothing to fix.** The
boundary jitter does not exist. It was verified against the two CARP runs that
motivated this request — the **same 13.76 Gbp genome**, same command
(`dante ... -c 32`), same DANTE 0.2.11, Viridiplantae_v4.0:

| Metric | Value |
| --- | --- |
| GFF3 records compared (`run-000097` vs `run-000117`) | 4,165,836 |
| Raw differing lines | 64,351 (1.5 %) |
| Lines differing in coordinates / score / strand (cols 1–8) | **0** |
| Lines with a genuinely different `Region_Hits` *set* | **0** |
| Differences remaining after applying 0.2.12's canonical sort | **0** |

Every one of the 64,351 differences is `Region_Hits_Classifications` order
(plus, on 3,339 `Ambiguous_domain` records, the `/`-joined `domain_type` in
`Name=`). Not one of 4.16 M records differs in start/end/score/strand, so there
is **no boundary jitter**. The first differing record is exactly the example in
this request — same locus `47747–47965`, identical `Best_Hit=…|SIRE`, only the
hit-list order swaps (`Lyco[195bp]` vs `SIRE[207bp]` as the leading entry); the
apparent "top hit change" is just the first element of a reordered list, not a
real classification change.

This matches the static analysis: `lastal` is deterministic (byte-identical
across runs, even with `-P` threads), and no code path routes hash/thread order
into a coordinate. Applying exactly what 0.2.12 does (sort `Region_Hits` + sort
`domain_type`) collapses all 64,351 differences to **0** — i.e. 0.2.12 makes the
two production 0.2.11 runs byte-identical. Fix #1 is therefore the complete
resolution; there is no separate boundary-determinism work to do.
