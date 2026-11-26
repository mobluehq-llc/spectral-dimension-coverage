#!/bin/bash
# BD chain on deep (WSL). Resumable: every mode skips done keys. Never kills anything. Usage: bash bd_chain.sh NCPU
. ~/hp/p22/env.sh; cd ~/hp/BD || exit 1; mkdir -p logs results
N=${1:-10}
for M in harness a b c d e f; do
  CMD="nice -n 10 python3 rt_cpbox.py $M --ncpu $N"
  echo "$(date -u +%FT%TZ) $CMD" >> logs/cmdline.log
  $CMD > logs/$M.log 2>&1 < /dev/null
  rc=$?; echo "$(date -u +%FT%TZ) mode $M exit $rc" >> logs/chain.log
  if [ $M = harness ] && ! grep -q '"pass": true' results/bd_harness_summary.json 2>/dev/null; then echo "$(date -u +%FT%TZ) HARNESS GATE FAILED, chain stopped" >> logs/chain.log; exit 3; fi
done
echo "$(date -u +%FT%TZ) chain complete" >> logs/chain.log
