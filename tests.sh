#!/bin/bash
# DANTE functional + determinism test.
#
#   1. Determinism gate: two runs of the same input with different
#      PYTHONHASHSEED must produce byte-identical domain output. This
#      guards the deterministic-output fix (canonical Region_Hits_
#      Classifications ordering + deterministic best-hit tie-break).
#   2. Correctness: the output body must match the committed reference
#      fixture (## comment/header lines are ignored, since they carry a
#      volatile git hash and version string).
#
# Runs against the source tree by default. Set DANTE to test an installed
# executable instead, e.g. `DANTE=dante ./tests.sh` for a conda-installed
# package. The env must have DANTE's runtime deps (numpy, biopython) and
# `lastal` on PATH.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

DANTE="${DANTE:-$ROOT/dante}"
DB="${DANTE_TEST_DB:-Viridiplantae_v4.0}"
INPUT="$ROOT/test-data/GEPY_test_long_1.fa"
REF="$ROOT/test-data/GEPY_test_long_1_output_unfiltered.gff3"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

strip_header() { grep -v '^##' "$1" || true; }

echo "== DANTE determinism + fixture test =="
echo "dante:    $DANTE"
echo "database: $DB"
echo "input:    $INPUT"

# 1) Determinism: two runs with different hash seeds must be byte-identical.
PYTHONHASHSEED=1   "$DANTE" -q "$INPUT" -o "$WORK/run1.gff3" -D "$DB" -c 1 >"$WORK/run1.log" 2>&1
PYTHONHASHSEED=999 "$DANTE" -q "$INPUT" -o "$WORK/run2.gff3" -D "$DB" -c 1 >"$WORK/run2.log" 2>&1
strip_header "$WORK/run1.gff3" > "$WORK/run1.body"
strip_header "$WORK/run2.gff3" > "$WORK/run2.body"
if ! cmp -s "$WORK/run1.body" "$WORK/run2.body"; then
  echo "FAIL: output is not deterministic across PYTHONHASHSEED"
  diff "$WORK/run1.body" "$WORK/run2.body" | head -20
  exit 1
fi
echo "PASS: deterministic across hash seeds ($(wc -l < "$WORK/run1.body") records)"

# 2) Correctness: output body matches the committed reference (## ignored).
strip_header "$REF" > "$WORK/ref.body"
if ! cmp -s "$WORK/run1.body" "$WORK/ref.body"; then
  echo "FAIL: output differs from reference $REF"
  echo "--- diff (< reference vs > current) ---"
  diff "$WORK/ref.body" "$WORK/run1.body" | head -40
  exit 1
fi
echo "PASS: output matches reference fixture"
echo "ALL TESTS PASSED"
