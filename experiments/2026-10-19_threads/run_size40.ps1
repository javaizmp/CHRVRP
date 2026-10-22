$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$root = Split-Path -Parent (Split-Path -Parent $here)
$py   = "$root\.venv\Scripts\python.exe"
$out  = "$here\bp_size40"
$work = "$out\work"              # bp_chrvrp.py writes <instance>.bp.json next to the instance, so run on copies
New-Item -ItemType Directory -Force $work | Out-Null
# warm the numba cache so compile time is not charged to the first instance
Copy-Item "$root\data\instances\small\instance20S2V4Y6PxLxn20.tsp" $work
& $py "$root\src\bp_chrvrp.py" "$work\instance20S2V4Y6PxLxn20.tsp" 6 --quiet | Out-Null
Remove-Item -Force "$work\instance20S2V4Y6PxLxn20.tsp*"
foreach ($t in 6, 1) {
    foreach ($n in 0..19) {
        $inst = "$work\instance_size_40_inst_$n.tsp"
        Copy-Item "$root\data\instances\size40\instance_size_40_inst_$n.tsp" $work
        $line = & $py "$root\src\bp_chrvrp.py" $inst 8 --quiet --threads $t
        "t=$t $line" | Out-File -Append -Encoding utf8 "$out\log.txt"
        Move-Item -Force "$inst.bp.json" "$out\inst${n}_t$t.json"
        Remove-Item -Force $inst
    }
}
Remove-Item -Force $work
'DONE' | Out-File -Append -Encoding utf8 "$out\log.txt"
