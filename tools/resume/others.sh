#!/bin/bash
# test_wsl_kill.ps1's refusal check, run inside each distro it would stop: prints every user process other
# than the resume dummy's session (its job, lock holder and their children), the user systemd's and this
# check's own. Output of only the CHECKED line means stopping the distro kills nothing but the dummy.
d=$(pgrep -f 'dummy_job[.]sh' | head -1)
ds=""
[ -n "$d" ] && ds=$(ps -o sid= -p "$d" | tr -d ' ')
me=$(ps -o sid= -p $$ | tr -d ' ')
ps -eo user=,sid=,etime=,args= | awk -v ds="$ds" -v me="$me" '
  $1 == "root" || $1 ~ /^(systemd|message|syslog|polkitd)/ { next }
  $2 == ds || $2 == me { next }
  / \(sd-pam\)$/ || /systemd --user/ { next }
  { print }'
# The caller goes ahead only on this line, so a check that fails to run refuses.
echo "CHECKED $(hostname)"
