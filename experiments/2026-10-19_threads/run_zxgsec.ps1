$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$root = Split-Path -Parent (Split-Path -Parent $here)
$py   = "$root\.venv\Scripts\python.exe"
$out  = "$here\zxgsec"
New-Item -ItemType Directory -Force $out | Out-Null
$env:CHRVRP_TIMELIMIT = 1800
foreach ($name in 'instance20S2V4Y6PxLxn20', 'instance20S2V4Y6PxLxn35') {
    foreach ($t in 1, 6) {
        $dir = "$out\${name}_t$t"
        New-Item -ItemType Directory -Force $dir | Out-Null
        Copy-Item "$root\data\instances\small\$name.tsp" $dir
        Set-Location $dir
        $env:CHRVRP_THREADS = $t
        & $py "$root\src\expt1_fast.py" "$name.tsp" zxgsec 6 *> stdout.txt
        "$name t=$t exit=$LASTEXITCODE" | Out-File -Append -Encoding utf8 "$out\progress.txt"
    }
}
'DONE' | Out-File -Append -Encoding utf8 "$out\progress.txt"
