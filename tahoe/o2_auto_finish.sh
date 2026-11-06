#!/bin/bash
# Wait for the Orion stream to end, then run the HEK293T final step.
cd "$(dirname "$0")"
until grep -q "all done" orion_o2_HEK293T.log 2>/dev/null; do
  pgrep -f run_o2.sh >/dev/null || { grep -q "finished HEK293T" orion_o2_HEK293T.log && break; echo "runner stopped before finishing" >&2; exit 1; }
  sleep 30
done
python3 o2_finalize.py HEK293T > o2_finalize_HEK293T.log 2>&1
grep -E "control cache|reference cells|median share|Spearman|^\(|^ *targets|^chronos" o2_finalize_HEK293T.log
