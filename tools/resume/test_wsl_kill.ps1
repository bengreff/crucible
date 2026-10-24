# The resume path's end-to-end test (tools/resume/README.md): stop WSL under the running dummy job, start the
# CrucibleResume task as boot would, and check the dummy resumes from its checkpoint. -Mode terminate stops
# the Ubuntu distro only; -Mode shutdown is wsl --shutdown, which stops every distro and the VM. It refuses
# if any user process other than the dummy runs in a distro it would stop, since stopping it would kill that job.
# -DryRun runs only the refusal check.
param([ValidateSet('terminate', 'shutdown')][string]$Mode = 'terminate', [switch]$DryRun)
$env:WSL_UTF8 = '1'
$D = '/home/greff/resume/dummy'
$distros = if ($Mode -eq 'shutdown') { @('Ubuntu', 'Ubuntu-22.04') } else { @('Ubuntu') }
foreach ($d in $distros) {
  # others.sh sits beside this script, which both distros see under /mnt/c.
  $out = @(wsl.exe -d $d -- bash /mnt/c/Users/Ben/crucible_resume/others.sh 2>&1 | Where-Object { "$_" -match '\S' -and "$_" -notmatch 'screen size is bogus' })
  $others = @($out | Where-Object { "$_" -notmatch '^CHECKED ' })
  if ($others.Count -gt 0 -or $out.Count -eq $others.Count) { "REFUSED: other processes run in $d; stopping it would kill them:"; $others; exit 2 }
}
if ($DryRun) { "dry run: nothing else runs in $($distros -join ', '); $Mode would go ahead"; exit 0 }
$before = (wsl.exe -d Ubuntu -- cat "$D/ckpt").Trim()
"$(Get-Date -Format 'HH:mm:ss') dummy's checkpoint before the stop: $before; mode $Mode"
if ($Mode -eq 'shutdown') { wsl.exe --shutdown } else { wsl.exe --terminate Ubuntu }
$t0 = Get-Date
do { Start-Sleep -Seconds 2; $state = (wsl.exe -l -v | Select-String '\sUbuntu\s').ToString() } until ($state -match 'Stopped' -or ((Get-Date) - $t0).TotalSeconds -gt 60)
"$(Get-Date -Format 'HH:mm:ss') after the stop: $($state.Trim())"
if ($state -notmatch 'Stopped') { 'FAILED: Ubuntu did not stop'; exit 1 }
Start-ScheduledTask -TaskName CrucibleResume
Start-Sleep -Seconds 30
Get-ScheduledTaskInfo -TaskName CrucibleResume | Format-List LastRunTime, LastTaskResult
$start = (wsl.exe -d Ubuntu -- tail -1 "$D/resume.log").Trim()
$a = [int](wsl.exe -d Ubuntu -- cat "$D/ckpt").Trim(); Start-Sleep -Seconds 10; $b = [int](wsl.exe -d Ubuntu -- cat "$D/ckpt").Trim()
"dummy's last start: $start"
"checkpoint 10 s apart: $a then $b"
if ($start -match 'start at (\d+)' -and [int]$Matches[1] -ge [int]$before - 1 -and $b -gt $a) { 'PASSED: resumed from its checkpoint and counting' } else { 'FAILED'; exit 1 }
