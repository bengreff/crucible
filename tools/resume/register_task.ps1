# Registers the Windows task CrucibleResume on backhouse (tools/resume/README.md): one minute after boot it
# starts WSL's Ubuntu as greff and runs resume.sh, which relaunches the registered jobs. It runs as Ben with
# logon type S4U, so it needs no logon and stores no password; it changes no account, logon or security setting.
# Run in an elevated PowerShell (an ssh session of an administrator is one). Re-running replaces the task.
$action = New-ScheduledTaskAction -Execute 'C:\Windows\System32\wsl.exe' `
  -Argument '-d Ubuntu -u greff --exec /bin/bash /home/greff/resume/resume.sh'
$trigger = New-ScheduledTaskTrigger -AtStartup
$trigger.Delay = 'PT1M'
$principal = New-ScheduledTaskPrincipal -UserId "$env:COMPUTERNAME\Ben" -LogonType S4U -RunLevel Limited
$settings = New-ScheduledTaskSettingsSet -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries `
  -ExecutionTimeLimit (New-TimeSpan -Minutes 10) -RestartCount 3 -RestartInterval (New-TimeSpan -Minutes 5) `
  -MultipleInstances IgnoreNew
Register-ScheduledTask -TaskName 'CrucibleResume' -Action $action -Trigger $trigger -Principal $principal `
  -Settings $settings -Description 'Crucible: relaunch registered WSL jobs from their checkpoints after a reboot (tools/resume/README.md)' -Force |
  Format-List TaskName,State
