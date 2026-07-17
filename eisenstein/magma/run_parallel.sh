#!/bin/bash
#
# Parallel regeneration of the Eisenstein CMF data tables (Option A: full
# DoEverythingNk2UpTo, split across many cores). Intended for a many-core box
# run it on any many-core server, inside tmux/screen (it is a long job).
#
# It runs three phases:
#   1. NJOBS parallel Magma workers, each doing DecomposeSpace on its slice
#      -> writes  <DATADIR>/Nk2_<B>.m.txt_<i>
#   2. concatenate + sort the slices by (N, k, char_orbit)  -> <DATADIR>/Nk2_<B>.m.txt
#      (sorting reproduces a single-process row order; data is identical either way)
#   3. one Magma process: dimension table + Format* -> the mf_*_<B>.m.txt tables
#
# Usage:
#   ./run_parallel.sh [B] [NJOBS] [DATADIR]
#     B        Nk^2 bound            (default 4000)
#     NJOBS    parallel workers      (default 64)
#     DATADIR  output folder, relative to eisenstein/magma (default ../data)
#
# Environment overrides:
#   MAGMA        magma binary            (default: magma on PATH)
#   DO_CONREY    true to (re)generate the conrey-labels file in phase 3
#                (default: auto -> true only if the file is missing)
#   KEEP_SLICES  true to keep the per-worker Nk2_<B>.m.txt_<i> files (default: false)
#   SORT_TMP     scratch dir for sort   (default: DATADIR)
#
set -euo pipefail

B="${1:-4000}"
NJOBS="${2:-64}"
DATADIR="${3:-../data}"
MAGMA="${MAGMA:-magma}"
KEEP_SLICES="${KEEP_SLICES:-false}"

cd "$(dirname "$0")"                 # -> eisenstein/magma
mkdir -p logs
INFILE="$DATADIR/Nk2_$B.m.txt"
CONREY="$DATADIR/Nk2_conrey_$B.m.txt"
SORT_TMP="${SORT_TMP:-$DATADIR}"

# Decide whether phase 3 must (re)generate the conrey-labels file.
if [ -z "${DO_CONREY:-}" ]; then
  if [ -f "$CONREY" ]; then DO_CONREY=false; else DO_CONREY=true; fi
fi

echo "=================================================================="
echo " Eisenstein regen: B=$B  NJOBS=$NJOBS  DATADIR=$DATADIR"
echo " MAGMA=$MAGMA  DO_CONREY=$DO_CONREY  (conrey file: $CONREY)"
echo "=================================================================="

# ---- Phase 1: parallel DecomposeSpace workers --------------------------------
echo "[$(date '+%F %T')] Phase 1: launching $NJOBS DecomposeSpace workers ..."
pids=()
for i in $(seq 0 $((NJOBS-1))); do
  "$MAGMA" B:=$B Jobs:=$NJOBS JobId:=$i DecomposeOnly:=true folder:="$DATADIR/" \
           run_eisenstein.m > "logs/worker_$i.log" 2>&1 &
  pids+=($!)
done
fail=0
for i in "${!pids[@]}"; do
  if ! wait "${pids[$i]}"; then
    echo "  !! worker $i FAILED (exit) -- see logs/worker_$i.log"; fail=1
  fi
done
# Magma reports internal errors without a non-zero exit; scan logs too.
if grep -lE "^(Runtime error|Magma:|>> )|error|Error|Killed" logs/worker_*.log >/dev/null 2>&1; then
  echo "  !! possible error markers found in worker logs:"; grep -lE "error|Error|Killed|Runtime" logs/worker_*.log || true
  fail=1
fi
[ $fail -eq 0 ] || { echo "Aborting after phase 1 (worker failure)."; exit 1; }
echo "[$(date '+%F %T')] Phase 1 done."

# ---- Phase 2: assemble the infile from the worker slices ---------------------
echo "[$(date '+%F %T')] Phase 2: concatenate + sort slices -> $INFILE"
slices=()
for i in $(seq 0 $((NJOBS-1))); do
  s="${INFILE}_$i"
  [ -f "$s" ] || { echo "  !! missing slice $s"; exit 1; }
  slices+=("$s")
done
cat "${slices[@]}" | LC_ALL=C sort -t: -k1,1n -k2,2n -k3,3n -T "$SORT_TMP" > "$INFILE"
echo "  assembled $(wc -l < "$INFILE") space-lines into $INFILE"
if [ "$KEEP_SLICES" != "true" ]; then rm -f "${slices[@]}"; fi

# ---- Phase 3: dimension table + Format* (single process) ---------------------
echo "[$(date '+%F %T')] Phase 3: dimension table + Format* ..."
"$MAGMA" B:=$B FormatOnly:=true folder:="$DATADIR/" do_conrey:=$DO_CONREY \
         run_eisenstein.m > logs/format.log 2>&1 || {
  echo "  !! Format phase failed -- see logs/format.log"; exit 1; }
echo "[$(date '+%F %T')] Phase 3 done."

# ---- Phase 4: manifest -------------------------------------------------------
MANIFEST="$DATADIR/MANIFEST_Nk2_$B.txt"
{
  echo "Eisenstein CMF data regeneration"
  echo "generated-by : eisenstein/magma/$(basename "$0") (B=$B, NJOBS=$NJOBS)"
  echo "host         : $(hostname)"
  echo "date         : $(date '+%F %T %Z')"
  echo "git-commit   : $(git -C .. rev-parse HEAD 2>/dev/null || echo '?')"
  echo ""
  echo "Tables (load these into the mf_*_eis Postgres tables):"
  for f in "$DATADIR"/Nk2_mf_*_"$B".m.txt "$DATADIR"/Nk2_dim"$B".m.txt; do
    [ -f "$f" ] && printf "  %-48s %12s bytes  %8s lines\n" \
        "$(basename "$f")" "$(wc -c < "$f")" "$(wc -l < "$f")"
  done
} | tee "$MANIFEST"
echo "[$(date '+%F %T')] All done. Manifest: $MANIFEST"
