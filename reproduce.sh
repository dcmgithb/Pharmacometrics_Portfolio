#!/usr/bin/env bash
# reproduce.sh — rebuild every analysis report from scratch.
#
# Prerequisites: R 4.4.1, Python >= 3.11, uv, Quarto >= 1.5
# Usage:  bash reproduce.sh
# Option: BOOTSTRAP_N_BOOT=50 bash reproduce.sh   (faster smoke-test run)
#
# The script is intentionally simple: restore envs → run scripts → render.
# Each project's scripts must be idempotent (overwrite their own outputs).

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ── Step 1: Restore R package library ────────────────────────────────────────
# renv::restore() reads renv.lock and installs pinned package versions into
# renv/library/, isolated from the system library.
log "Restoring R packages from renv.lock ..."
Rscript -e "if (!requireNamespace('renv', quietly=TRUE)) install.packages('renv', repos='https://cloud.r-project.org'); renv::restore(prompt=FALSE)"
log "R packages ready."

# ── Step 2: Sync Python virtual environment ───────────────────────────────────
# uv sync reads pyproject.toml and creates .venv/ with pinned versions.
log "Syncing Python dependencies via uv ..."
uv sync
log "Python packages ready."

# ── Step 3: Run R-based projects (01, 02, 03, 05) ────────────────────────────
R_PROJECTS=(
  "projects/01_popPK_nlmixr2"
  "projects/02_PKPD_emax"
  "projects/03_preclinical_to_human"
  "projects/05_trial_simulation"
)

for proj in "${R_PROJECTS[@]}"; do
  log "=== $proj ==="
  (
    cd "$REPO_ROOT/$proj"
    Rscript 01_simulate_data.R
    Rscript 02_fit_model.R
    Rscript 03_diagnostics.R
  )
  log "=== $proj done ==="
done

# ── Step 4: Run Python-based project (04) ────────────────────────────────────
log "=== projects/04_PBPK_allometry (Python) ==="
(
  cd "$REPO_ROOT/projects/04_PBPK_allometry"
  uv run python 01_simulate_data.py
  uv run python 02_fit_model.py
  uv run python 03_diagnostics.py
)
log "=== projects/04_PBPK_allometry done ==="

# ── Step 5: Render Quarto site ────────────────────────────────────────────────
# freeze: auto in _quarto.yml means only changed .qmd files re-execute.
# On a first run, all files execute. On subsequent runs, only changed ones do.
log "Rendering Quarto site to docs/ ..."
quarto render
log "Site rendered. Open docs/index.html to preview locally."
