#!/bin/bash
# Phase O2 runner: HCT116 then HEK293T; each line resumes from its checkpoint, and a crashed run is restarted.
cd "$(dirname "$0")"
for line in HCT116 HEK293T; do
  for attempt in $(seq 1 50); do
    python3 o2_orion_stream.py $line >> orion_o2_$line.log 2>&1 && break
    echo "$(date +%T) restart $line after failure (attempt $attempt)" >> orion_o2_$line.log; sleep 60
  done
done
echo "$(date +%T) all done" >> orion_o2_HEK293T.log
