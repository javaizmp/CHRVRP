# CHRVRP project – notes for Claude

Code, instances, experiments and manuscript for the CHR-VRP paper. Layout, setup and commands are in README.md.

## Environment
- Windows only. Use `.venv\Scripts\python.exe` (Python 3.10). Do not use the system `python` (3.14, no CPLEX).
- CPLEX is the full Studio 22.1.1 API, installed via its `setup.py`, never `pip install cplex` (that is the size-limited Community Edition). See README "Setup".
- No WSL, no git, no Perl. Build the paper with `pdflatex main` / `bibtex main` / `pdflatex` ×2, not latexmk.

## Paths
- Every script uses paths relative to the project root. Never hardcode `D:\...` or `C:\Users\...` paths.
- `jobs.csv` stores instance paths relative to the root; `run_jobs.ps1` resolves them.

## Running solvers
- Solvers write output to their working directory:
  `<inst>.tsprs` (result summary, parsed into results.csv), `<inst>.tsplog` (CPLEX log),
  `agg<name>.rs` / `agg_*log` (expt_agg.py), `<inst>.bp.json` (bp_chrvrp.py, written next to the instance file).
- So run them via `experiments\2026-10\run_jobs.ps1` (each job gets `runs\<id>\`) or the scripts in
  `experiments\2026-10-19_threads\`. For one-off tests, copy the instance to a scratch/temp folder and run there.
  Never run a solver inside `src\` or `data\`.
- To re-check results without touching `runs\`: `run_jobs.ps1 -Only '<pattern>' -RunsDir <scratch>\runs`, then compare with `results.csv`.
- `results.csv` has exactly one row per job (the latest run), matching `runs\<id>\result.json`. To re-run a job, delete or move
  its `runs\<id>\` folder; the runner then replaces that job's row. Older rows from the 2026-10 re-runs are in `superseded\results_superseded.csv`.
- Output files inside `experiments\` are raw results. Keep them. Elsewhere they are strays (.gitignore reflects this).

## Keeping the project clean
- No versioned copies (`script_22May.py`, `main_backup_<date>.tex`): edit in place. Use git for history once it is installed.
- One copy of each instance under `data\instances\`. Don't duplicate instance sets.
- New experiment campaigns go in `experiments\<yyyy-mm>[_topic]\` with their own make/run scripts.
- Paper: `paper\main.tex` is the manuscript (formerly `main_corrections_bp.tex`). Figure PDFs sit next to `main.tex`. LaTeX build files are not kept.
- **Never edit `paper\main.tex` above the line `%%%%%%%JAVAIZ EDITED%%%%%%%`.** That part, including the preamble, has been finalised by the user.
  Make all edits below the marker. If a change seems to need something above it (a package, a macro, a fix in the intro or abstract),
  ask the user rather than editing it. Before editing, find where the marker is now; the user moves it as more of the paper is finalised.
- New text in the paper uses `\textcolor{revpink}{...}` (tables: `\color{revpink}` inside the `table` environment), matching the branch-and-price revisions.

## GitHub
- Repository is public on GitHub (account javaizmp). Commit as the GitHub noreply identity (set in this repo's git config).
- `paper/` is gitignored: the manuscript stays local until accepted. Never force-add it or paste its contents into committed files.

## Copying the project
- Copying the folder works anywhere (relative paths). `.venv` is the exception: it works when copied on this PC
  (use `.venv\Scripts\python.exe -m pip`, the pip/activate launchers keep the old path), but not on another PC.
  There, leave `.venv` out and rebuild it (README "Setup"), or prefer `git clone` and then rebuild `.venv`.
- `paper/` is not on GitHub. A folder copy is its only backup, so remind the user to include it.
- Copies share the GitHub remote. When there are several working copies, push/pull before working so they don't diverge.
- To check a copy: `run_jobs.ps1 -Only 'P1_S1_bp' -RunsDir <temp>\runs` → optimal, obj 117.89.

## History
- Migrated on 2026-10-22 from `D:\PhD\Courses\RAM_CHRBEND\Python Implementation`. KPI3 scripts, the 2023 Benders
  experiments (`Final Experiments` etc.), old .tex versions and backups were left there (to be archived) and are not part of this project.
- Migration check: 32 jobs (all four solvers, S1/S2) reproduced results.csv exactly, and the paper compiled identically.
