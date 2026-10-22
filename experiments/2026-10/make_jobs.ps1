# Builds jobs.csv for the 2026-10 experiment campaign. Re-running overwrites jobs.csv
# but never touches results (run_jobs.ps1 skips jobs whose result file exists).
#
# Common settings: 6 threads, 3600 s limit, all lower bounds (CHRVRP_LB=1234) on
# every arc-based model, vehicles optional (paper da5 "<= 1") everywhere.
#
# Instance paths are written relative to the project root so jobs.csv stays portable.
$here = Split-Path -Parent $MyInvocation.MyCommand.Path

# S1-S20 -> instance file (mapped by the optimum in the paper's Table 1)
$small = [ordered]@{
    S1='instance10S2V3Y4PxLxn20'; S2='instance10S2V3Y4PxLxn5';  S3='instance10S2V3Y4PxLxn22'
    S4='instance10S2V3Y4PxLxn23'; S5='instance10S2V3Y4PxLxn24'; S6='instance10S2V3Y4PxLxn25'
    S7='instance10S2V3Y4PxLxn26'; S8='instance10S2V3Y4PxLxn27'; S9='instance10S2V3Y4PxLxn10'
    S10='instance10S2V3Y4PxLxn29'
    S11='instance20S2V4Y6PxLxn20'; S12='instance20S2V4Y6PxLxn35'; S13='instance20S2V4Y6PxLxn22'
    S14='instance20S2V4Y6PxLxn23'; S15='instance20S2V4Y6PxLxn24'; S16='instance20S2V4Y6PxLxn25'
    S17='instance20S2V4Y6PxLxn26'; S18='instance20S2V4Y6PxLxn27'; S19='instance20S2V4Y6PxLxn50'
    S20='instance20S2V4Y6PxLxn29'
}
$jobs = New-Object System.Collections.Generic.List[object]
function Add-Job($phase, $label, $inst, $q, $method, $lb = '', $beta = 0, $lprel = 0) {
    $id = '{0}_{1}_{2}' -f $phase, $label, $method
    if ($lprel) { $id += '_lp' + $(if ($lb) { $lb } else { 'none' }) }
    $jobs.Add([pscustomobject]@{ id = $id; phase = $phase; label = $label; instance = $inst; Q = $q
        method = $method; LB = $lb; BETATIGHT = $beta; LPREL = $lprel })
}
$smallJobs = foreach ($k in $small.Keys) { [pscustomobject]@{ label = $k; inst = "data\instances\small\$($small[$k]).tsp"; q = $(if ([int]$k.Substring(1) -le 10) { 4 } else { 6 }) } }
$s30 = 0..19 | ForEach-Object { [pscustomobject]@{ label = "d30_$_"; inst = "data\instances\size30\V5Q7\instance_size_30_inst_$_.tsp"; q = 7 } }
$s40 = 0..19 | ForEach-Object { [pscustomobject]@{ label = "d40_$_"; inst = "data\instances\size40\instance_size_40_inst_$_.tsp"; q = 8 } }

# Block A (minutes): B&P everywhere, and LP relaxations for "Best LP Rel" and the LP/LB figures
foreach ($s in $smallJobs + $s30 + $s40) { Add-Job ($(if ($s.label -like 'S*') { 'P1' } else { 'P2' })) $s.label $s.inst $s.q 'bp' }
foreach ($s in $smallJobs) {
    foreach ($lb in '', '1234') {
        Add-Job 'P1' $s.label $s.inst $s.q 'zxdefault' $lb 0 1
        Add-Job 'P1' $s.label $s.inst $s.q 'azf' $lb 0 1
    }
    foreach ($lb in '', '1', '2', '3', '4', '1234') { Add-Job 'P1' $s.label $s.inst $s.q 'zfdefault' $lb 0 1 }
}
# Block B (hours): fast arc models on S1-S20
foreach ($s in $smallJobs) { Add-Job 'P1' $s.label $s.inst $s.q 'zxgsec' '1234' }
foreach ($s in $smallJobs) { Add-Job 'P1' $s.label $s.inst $s.q 'zxdefault' '1234' }
# Block C (~35 h): best arc model vs B&P on Size 30 (V5Q7) and Size 40
foreach ($s in $s30 + $s40) { Add-Job 'P2' $s.label $s.inst $s.q 'zxgsec' '1234' }
# Block D (~50 h): slow arc models on S1-S20 (S11-S20 expected to hit the limit)
foreach ($s in $smallJobs) { Add-Job 'P1' $s.label $s.inst $s.q 'zftbend' '1234' 1 }
foreach ($s in $smallJobs) { Add-Job 'P1' $s.label $s.inst $s.q 'zftbend_nobeta' '1234' 0 }
foreach ($s in $smallJobs) { Add-Job 'P1' $s.label $s.inst $s.q 'zfdefault' '1234' }
foreach ($s in $smallJobs) { Add-Job 'P1' $s.label $s.inst $s.q 'azf' '1234' }

$jobs | Export-Csv -NoTypeInformation -Encoding UTF8 "$here\jobs.csv"
"{0} jobs written to {1}\jobs.csv" -f $jobs.Count, $here
