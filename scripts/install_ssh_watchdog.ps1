<#
.SYNOPSIS
Register (or remove) a per-user scheduled task that probes the cluster's SSH
banner every few minutes.

.DESCRIPTION
The task runs scripts/ssh_watchdog.py with the host supplied at install time and
appends one JSON line per probe to %LOCALAPPDATA%\rff\ssh_watchdog.log. It is
user-level (no administrator rights) and never restarts the VPN client or edits
routes; it exists so a silent tunnel loss is visible in the log instead of
looking like an idle server.

.EXAMPLE
powershell -File scripts/install_ssh_watchdog.ps1 -SshHost 10.0.0.5 -EveryMinutes 5

.EXAMPLE
powershell -File scripts/install_ssh_watchdog.ps1 -Uninstall
#>
param(
    [string]$SshHost,
    [int]$Port = 22,
    [int]$EveryMinutes = 5,
    [string]$PythonPath,
    [switch]$Uninstall
)

$ErrorActionPreference = "Stop"
$TaskName = "RFF SSH Watchdog"
$Root = Split-Path -Parent $PSScriptRoot
$Script = Join-Path $Root "scripts\ssh_watchdog.py"

if ($Uninstall) {
    if (Get-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue) {
        Unregister-ScheduledTask -TaskName $TaskName -Confirm:$false
        Write-Output "removed scheduled task '$TaskName'"
    } else {
        Write-Output "no scheduled task '$TaskName' to remove"
    }
    return
}

if (-not $SshHost) { throw "-SshHost is required unless -Uninstall is given" }
if (-not (Test-Path $Script)) { throw "watchdog script not found: $Script" }

# The Windows Store alias `python.exe` is a 0-byte stub on many machines, so a
# resolved path is not enough: each candidate must actually run.
function Test-Python([string]$Exe, [string[]]$Prefix) {
    try {
        $out = & $Exe @Prefix --version 2>&1
        return ($LASTEXITCODE -eq 0 -and "$out" -match "Python 3")
    } catch { return $false }
}

$Python = $null
$PythonPrefix = @()
if ($PythonPath) {
    if (-not (Test-Python $PythonPath @())) { throw "-PythonPath does not run: $PythonPath" }
    $Python = $PythonPath
} else {
    $Py = (Get-Command py -ErrorAction SilentlyContinue).Source
    if ($Py -and (Test-Python $Py @("-3"))) { $Python = $Py; $PythonPrefix = @("-3") }
    if (-not $Python) {
        $Py3 = (Get-Command python3 -ErrorAction SilentlyContinue).Source
        if ($Py3 -and (Test-Python $Py3 @())) { $Python = $Py3 }
    }
    if (-not $Python) { throw "no working Python 3 found; pass -PythonPath" }
}

$ArgPrefix = if ($PythonPrefix.Count) { ($PythonPrefix -join " ") + " " } else { "" }
$Action = New-ScheduledTaskAction -Execute $Python `
    -Argument "$ArgPrefix`"$Script`" --host $SshHost --port $Port --quiet" `
    -WorkingDirectory $Root

$Trigger = New-ScheduledTaskTrigger -Once -At (Get-Date).AddMinutes(1) `
    -RepetitionInterval (New-TimeSpan -Minutes $EveryMinutes)

$Settings = New-ScheduledTaskSettingsSet -StartWhenAvailable `
    -MultipleInstances IgnoreNew -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries

$Principal = New-ScheduledTaskPrincipal -UserId "$env:USERDOMAIN\$env:USERNAME" `
    -LogonType Interactive -RunLevel Limited

Register-ScheduledTask -TaskName $TaskName -Action $Action -Trigger $Trigger `
    -Settings $Settings -Principal $Principal -Force | Out-Null

Write-Output "registered '$TaskName': every $EveryMinutes min, host $SshHost port $Port"
Write-Output "log: $env:LOCALAPPDATA\rff\ssh_watchdog.log"
