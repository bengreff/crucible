#!/bin/bash
# Backhouse's resume-after-reboot runner (tools/resume/README.md). With no argument, it relaunches every job
# registered in $RESUME_HOME/jobs/*.job that is not done, not held and not running; the Windows task
# CrucibleResume runs it at boot. With a name, it launches that one job the same way, which is how a job is
# first started, so every registered job runs under its lock.
R=${RESUME_HOME:-/home/greff/resume}
mkdir -p "$R/jobs"
log() { echo "$(date '+%F %T %Z') $*" >> "$R/resume.log"; }
log "resume.sh ${1:-all}: up $(cut -d' ' -f1 /proc/uptime) s, load $(cut -d' ' -f1 /proc/loadavg)"
if [ -n "$1" ]; then set -- "$R/jobs/$1.job"; else set -- "$R"/jobs/*.job; fi
for f in "$@"; do
  [ -e "$f" ] || { log "no job file $f"; continue; }
  (
    NAME="" DIR="" RESUME="" DONE="" NICE=10
    . "$f"
    if [ -z "$NAME" ] || [ -z "$DIR" ] || [ -z "$RESUME" ] || [ -z "$DONE" ]; then log "$f: NAME, DIR, RESUME or DONE unset; skipped"; exit 0; fi
    if [ -e "$DONE" ]; then log "$NAME: done ($DONE); skipped"; exit 0; fi
    if [ -e "$DIR/.resume_hold" ]; then log "$NAME: held ($DIR/.resume_hold); skipped"; exit 0; fi
    cd "$DIR" || { log "$NAME: no directory $DIR; skipped"; exit 0; }
    # The job holds this lock for as long as it runs; the kernel drops it when the process dies, WSL's end included.
    if ! flock -n .resume.lock true; then log "$NAME: running (lock held); skipped"; exit 0; fi
    # A login shell, so the job sees what an ssh "wsl bash -l" launch sees (PATH, the venvs, Cantera's paths).
    setsid nohup flock -n .resume.lock nice -n "$NICE" bash -lc "$RESUME" >> resume.log 2>&1 < /dev/null &
    log "$NAME: launched in $DIR, pid $!: $RESUME"
  )
done
exit 0
