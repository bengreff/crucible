#!/bin/bash
# The resume path's test job (tools/resume/README.md): counts to $1 (default 600), one step a second. It
# checkpoints as the engine does (CRUCIBLE_CHECKPOINT, inquiry's interface): write ckpt.tmp, sync, hard-link
# the old ckpt to ckpt.prev, rename. On start it resumes from ckpt, or ckpt.prev if ckpt is unreadable, or
# from 0. Its done-file, done, is written only when it reaches the end.
N=${1:-600}
k=0
for c in ckpt ckpt.prev; do
  if [ -f "$c" ] && grep -qx '[0-9][0-9]*' "$c"; then k=$(cat "$c"); break; fi
done
echo "$(date '+%F %T %Z') start at $k of $N (pid $$)"
while [ "$k" -lt "$N" ]; do
  sleep 1
  k=$((k + 1))
  echo "$k" > ckpt.tmp && sync ckpt.tmp && { [ ! -f ckpt ] || ln -f ckpt ckpt.prev; } && mv ckpt.tmp ckpt
done
echo "$(date '+%F %T %Z') reached $N"
touch done
