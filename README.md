# CHRVRP

Code, instances, experiments and manuscript for the CHR-VRP paper
(Benders arc models and branch-and-price).

## Layout

```
src/                     solvers (run with the instance file in the working directory)
  bp_chrvrp.py           branch-and-price                  bp_chrvrp.py <inst> <Q> [--threads n --timelimit s]
  expt1_fast.py          disaggregated arc models          expt1_fast.py <inst> <mode> <Q>   (zxgsec, zxdefault, zftbend, zfdefault, ...)
  expt_agg.py            aggregated model (azf)            expt_agg.py <inst> 0 <Q>
  chrvrp_io.py           .tsp instance reader
data/
  instances/small/       S1-S20 (paper Table 1; mapping in experiments/2026-10/make_jobs.ps1)
  instances/size30/      Size 30, one folder per V/Q setting (experiments use V5Q7)
  instances/size40/      Size 40
  generator/             instance generator script
experiments/
  2026-10/               main campaign: make_jobs.ps1 -> jobs.csv -> run_jobs.ps1 -> runs/, results.csv
  2026-10/superseded/    earlier B&P runs replaced by runs/
  2026-10-19_threads/    1 vs 6 thread comparison (B&P size 40, zxgsec S11/S12)
paper/                   manuscript (local only; not in the public repository until the paper is accepted)
```

All scripts use paths relative to the project root, so the folder can be moved or copied as a whole.

## Setup (Windows)

```powershell
py -3.10 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

CPLEX: install IBM ILOG CPLEX Studio 22.1.1 (academic licence), then install its Python API into
the venv. Its `setup.py` looks for `..\cplex\python\3.10`, so run it from the Studio folder layout:

```powershell
cd "C:\Program Files\IBM\ILOG\CPLEX_Studio2211\python"
D:\CHRVRP\.venv\Scripts\python.exe setup.py install
```

(If Program Files is not writable, copy `CPLEX_Studio2211\python` and `CPLEX_Studio2211\cplex\python\3.10`
to a temp folder keeping the same relative layout, and run it there.)

## Copying or moving the project

All paths are relative, so the folder works from any location. Only `.venv` depends on the machine:

- **Same PC, new location:** copy the whole folder. `.venv\Scripts\python.exe` still works because it uses
  this PC's Python 3.10 and has CPLEX inside. The `pip.exe`/`activate` launchers still point to the old path,
  so use `.venv\Scripts\python.exe -m pip ...`.
- **Another PC, USB stick or backup:** copy everything except `.venv`, then rebuild it as in *Setup* above.
- **Another PC, preferred:** `git clone https://github.com/javaizmp/CHRVRP`, rebuild `.venv`, and copy `paper/`
  across separately (it is not on GitHub).

A copy includes `.git`, so it is linked to the same GitHub repository. If you work in more than one copy,
`git push` from one and `git pull` in the other before starting, so they don't diverge.

Check that a copy works:

```powershell
cd experiments\2026-10
.\run_jobs.ps1 -Only 'P1_S1_bp' -RunsDir $env:TEMP\check\runs      # should finish in seconds, obj 117.89
```

## Running experiments

```powershell
experiments\2026-10\make_jobs.ps1                      # rebuild jobs.csv
experiments\2026-10\run_jobs.ps1                       # run all jobs not yet done (resumable)
experiments\2026-10\run_jobs.ps1 -Only 'P1_S1_*'       # subset
experiments\2026-10\run_jobs.ps1 -Only 'P1_S1_*' -RunsDir C:\temp\check\runs   # re-run elsewhere without touching runs/
```

Create a file named `STOP` in `experiments\2026-10` to stop after the current job.

## Building the paper

From `paper/`: `pdflatex main`, `bibtex main`, `pdflatex main`, `pdflatex main`
(`latexmk` also works if Perl is installed).
