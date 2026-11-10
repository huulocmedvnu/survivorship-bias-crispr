#!/bin/bash
# Stream X-Atlas/Orion for HCT116 then HEK293T. Each line resumes from its checkpoint after a failure.
HERE="$(cd "$(dirname "$0")" && pwd)"; LOG="$HERE/../../data/orion/logs"; mkdir -p "$LOG"
for line in HCT116 HEK293T; do
  for attempt in $(seq 1 50); do
    python3 "$HERE/11_orion_stream.py" $line >> "$LOG/stream_$line.log" 2>&1 && break
    echo "$(date +%T) restart $line after failure (attempt $attempt)" >> "$LOG/stream_$line.log"; sleep 60
  done
done
