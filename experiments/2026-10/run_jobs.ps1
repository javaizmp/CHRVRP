# Runs the jobs in jobs.csv one at a time.
#   * resumable: a job with runs\<id>\result.json is skipped
#   * pause: create a file named STOP next to this script; the runner exits after
#     the current job (delete STOP and start again to resume)
#   * like-for-like check: every 15 s the CPU use of all other processes is
#     sampled; a job whose other-process load is high is flagged in results.csv
#   * keeps Windows awake while running
param([string]$Only = '',         # optional wildcard on job id, e.g. 'P1_*_bp'
      [int]$TimeLimit = 3600,       # lower only for smoke tests
      [string]$RunsDir = '')        # default: .\runs (results.csv/progress.txt sit next to it)
$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$root = Split-Path -Parent (Split-Path -Parent $here)     # project root (instance paths in jobs.csv are relative to it)
$src  = "$root\src"
$py   = "$root\.venv\Scripts\python.exe"
$THREADS = 6
$TIMELIMIT = $TimeLimit
$GRACE = 900                 # hard kill at TIMELIMIT + GRACE (model build, LP dumps)
$FLAG_MEAN = 15              # flag if other processes used > 15% of one core on average
$FLAG_MAX = 100              # ... or more than one full core in any sample

Add-Type -Namespace W -Name P -MemberDefinition '[DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint f);'
[W.P]::SetThreadExecutionState([uint32]2147483649) | Out-Null    # ES_CONTINUOUS | ES_SYSTEM_REQUIRED

$runs = $(if ($RunsDir) { $RunsDir } else { "$here\runs" }); New-Item -ItemType Directory -Force $runs | Out-Null
$resCsv = "$(Split-Path -Parent $runs)\results$(if ($RunsDir) { '_' + (Split-Path -Leaf $runs) }).csv"
$progress = "$(Split-Path -Parent $runs)\progress$(if ($RunsDir) { '_' + (Split-Path -Leaf $runs) }).txt"
$jobs = Import-Csv "$here\jobs.csv"
if ($Only) { $jobs = $jobs | Where-Object { $_.id -like $Only } }

function Read-Lines($f) { if (Test-Path $f) { [string[]](Get-Content $f) } else { $null } }   # plain strings, so result.json holds values not PS objects

function Sample-Others([int[]]$exclude) {
    $rows = Get-CimInstance Win32_PerfFormattedData_PerfProc_Process -Property Name, IDProcess, PercentProcessorTime
    $o = $rows | Where-Object { $_.Name -notin '_Total', 'Idle' -and $_.IDProcess -notin $exclude -and $_.PercentProcessorTime -gt 0 }
    [pscustomobject]@{ total = [double](($o | Measure-Object PercentProcessorTime -Sum).Sum)
                       top = @($o | Sort-Object PercentProcessorTime -Descending | Select-Object -First 3 | ForEach-Object { '{0}:{1}' -f ($_.Name -replace '#\d+$', ''), $_.PercentProcessorTime }) }
}

$n = 0
foreach ($j in $jobs) {
    $n++
    if (Test-Path "$here\STOP") { "STOP file found, exiting before $($j.id)" | Out-File -Append -Encoding utf8 $progress; break }
    $dir = "$runs\$($j.id)"
    if (Test-Path "$dir\result.json") { continue }
    New-Item -ItemType Directory -Force $dir | Out-Null
    Copy-Item (Join-Path $root $j.instance) $dir -Force
    $file = Split-Path $j.instance -Leaf
    $base = [IO.Path]::GetFileNameWithoutExtension($file)

    $env:CHRVRP_THREADS = $THREADS; $env:CHRVRP_TIMELIMIT = $TIMELIMIT
    $env:CHRVRP_LB = $j.LB; $env:CHRVRP_BETATIGHT = $j.BETATIGHT
    $env:CHRVRP_LPREL = $(if ($j.LPREL -eq '1') { '1' } else { '' })
    $env:CHRVRP_FORCED_VEHICLES = ''; $env:CHRVRP_WRITELP = ''
    switch -Wildcard ($j.method) {
        'bp'      { $args_ = @("`"$src\bp_chrvrp.py`"", $file, $j.Q, '--threads', $THREADS, '--timelimit', $TIMELIMIT, '--quiet') }
        'azf'     { $args_ = @("`"$src\expt_agg.py`"", $file, 0, $j.Q) }
        'zftbend*' { $args_ = @("`"$src\expt1_fast.py`"", $file, 'zftbend', $j.Q) }
        default   { $args_ = @("`"$src\expt1_fast.py`"", $file, $j.method, $j.Q) }
    }
    $samples = New-Object System.Collections.Generic.List[double]
    $tops = @{}
    $s = Sample-Others @($PID)          # baseline just before launch, so short jobs get a reading too
    $samples.Add($s.total)
    foreach ($t in $s.top) { $nm = $t.Split(':')[0]; $tops[$nm] = 1 + [int]$tops[$nm] }
    $start = Get-Date
    "[{0:yyyy-MM-dd HH:mm}] {1}/{2} {3}" -f $start, $n, $jobs.Count, $j.id | Out-File -Append -Encoding utf8 $progress
    $p = Start-Process -FilePath $py -ArgumentList $args_ -WorkingDirectory $dir -NoNewWindow -PassThru `
                       -RedirectStandardOutput "$dir\stdout.txt" -RedirectStandardError "$dir\stderr.txt"
    $exclude = @($PID, $p.Id)
    $killed = $false
    while (-not $p.WaitForExit(15000)) {
        $s = Sample-Others $exclude
        $samples.Add($s.total)
        foreach ($t in $s.top) { $nm = $t.Split(':')[0]; $tops[$nm] = 1 + [int]$tops[$nm] }
        if (((Get-Date) - $start).TotalSeconds -gt $TIMELIMIT + $GRACE) { Stop-Process -Id $p.Id -Force; $killed = $true; break }
    }
    $wall = ((Get-Date) - $start).TotalSeconds

    $r = [ordered]@{ id = $j.id; phase = $j.phase; label = $j.label; method = $j.method; LB = $j.LB
                     BETATIGHT = $j.BETATIGHT; LPREL = $j.LPREL; Q = $j.Q; threads = $THREADS
                     status = ''; obj = ''; time = ''; gap = ''; bound = ''; nodes = ''; lp = ''
                     wall = [math]::Round($wall, 1); exit = $p.ExitCode; killed = $killed }
    if ($j.LPREL -eq '1') {
        $l = Read-Lines "$dir\${file}lprel"
        if ($l) { $r.lp = $l[0]; $r.status = $l[1] }
    } elseif ($j.method -eq 'bp') {
        if (Test-Path "$dir\$file.bp.json") {
            $b = Get-Content "$dir\$file.bp.json" -Raw | ConvertFrom-Json
            $r.status = $b.status; $r.obj = $b.obj; $r.time = [math]::Round($b.time, 2); $r.bound = $b.lb
            $r.gap = [math]::Round(($b.obj - $b.lb) / $b.obj, 6); $r.nodes = $b.nodes; $r.lp = $b.root_lb
        }
    } else {
        $rs = if ($j.method -eq 'azf') { Read-Lines "$dir\agg$base.rs" } else { Read-Lines "$dir\${file}rs" }
        if ($rs) { $r.obj = $rs[1]; $r.status = $rs[2]; $r.time = $rs[3]; $r.gap = $rs[4]; $r.bound = $rs[5]; if ($rs.Count -gt 6) { $r.nodes = $rs[6] } }
    }
    $mean = if ($samples.Count) { ($samples | Measure-Object -Average).Average } else { 0 }
    $max = if ($samples.Count) { ($samples | Measure-Object -Maximum).Maximum } else { 0 }
    $r.other_cpu_mean = [math]::Round($mean, 1)
    $r.other_cpu_max = [math]::Round($max, 1)
    $r.other_top = (($tops.GetEnumerator() | Sort-Object Value -Descending | Select-Object -First 3 | ForEach-Object { $_.Key }) -join ' ')
    $r.flag = @($(if ($mean -gt $FLAG_MEAN -or $max -gt $FLAG_MAX) { 'OTHER_LOAD' }),
                $(if ($killed) { 'KILLED' }),
                $(if (-not $r.status) { 'NO_RESULT' })) -join ' '
    $r.flag = $r.flag.Trim()
    $r.start = $start.ToString('s')
    [pscustomobject]$r | ConvertTo-Json | Out-File -Encoding utf8 "$dir\result.json"
    # one row per job: a re-run replaces the job's earlier row
    $prev = @(if (Test-Path $resCsv) { Import-Csv $resCsv | Where-Object { $_.id -ne $j.id } })
    @($prev) + [pscustomobject]$r | Export-Csv -NoTypeInformation -Encoding UTF8 $resCsv
    "    -> status={0} obj={1} time={2} wall={3:N0}s other_cpu={4}/{5} {6}" -f $r.status, $r.obj, $r.time, $wall, $r.other_cpu_mean, $r.other_cpu_max, $r.flag |
        Out-File -Append -Encoding utf8 $progress
}
[W.P]::SetThreadExecutionState([uint32]2147483648) | Out-Null    # back to normal
"[{0:yyyy-MM-dd HH:mm}] runner finished" -f (Get-Date) | Out-File -Append -Encoding utf8 $progress
