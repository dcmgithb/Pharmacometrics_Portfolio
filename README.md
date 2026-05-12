# Pharmacometrics Portfolio

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Quarto](https://img.shields.io/badge/Built_with-Quarto-blue)](https://quarto.org)
[![R](https://img.shields.io/badge/R-4.4.1-276DC3)](https://www.r-project.org)
[![Python](https://img.shields.io/badge/Python-3.11-3776AB)](https://www.python.org)

Five reproducible pharmacometrics analyses demonstrating core competencies
across the drug development pipeline: population PK modelling, sequential
PK/PD, preclinical-to-human translation, minimal PBPK, and clinical trial
simulation. The toolchain is entirely open-source (nlmixr2, rxode2, mrgsolve,
scipy); no NONMEM licence is required. Each analysis is rendered as a
self-contained HTML report that a reviewer can read in under five minutes.

> **All datasets in this repository are entirely simulated.**  
> No patient data, proprietary data, or licensed software (e.g., NONMEM) is
> used at any point. True parameter values are documented at the top of every
> `01_simulate_data.*` script so that parameter recovery can be verified.

---

## Projects

| # | Title | Clinical Question | Key Output | Tools |
|---|-------|-------------------|------------|-------|
| [01](projects/01_popPK_nlmixr2/index.qmd) | Population PK — 2-cmt oral | Does flat 100 mg QD maintain target trough across BW 40–120 kg, or is weight-based dosing needed? | Covariate-adjusted dosing decision with simulated trough distributions | nlmixr2 (SAEM), rxode2 |
| [02](projects/02_PKPD_emax/index.qmd) | Sigmoid Emax PK/PD | Which exposure metric best predicts a continuous biomarker, and what is the EC₅₀? | Recommended Phase 2 dose range with 90 % CI | nlmixr2, ggplot2 |
| [03](projects/03_preclinical_to_human/index.qmd) | Preclinical → Human | What FIH starting dose is supported by rat and dog PK with allometric scaling? | FIH dose with MABEL/NOAEL safety margin | nlmixr2, allometry |
| [04](projects/04_PBPK_allometry/index.qmd) | Minimal PBPK | Does allometric clearance in a minimal PBPK model reproduce human PK across age groups? | Pediatric (6 yo, 25 kg) dose recommendation | scipy, matplotlib |
| [05](projects/05_trial_simulation/index.qmd) | Clinical Trial Simulation | What sample size achieves ≥ 80 % power to estimate EC₅₀ within ± 30 %? | Power curve and operating characteristics table | rxode2, future.apply |

---

## Reproduce in One Command

**Prerequisites:** R 4.4.1, Python ≥ 3.11, [uv](https://docs.astral.sh/uv/),
[Quarto](https://quarto.org) ≥ 1.5.

```bash
git clone https://github.com/dcmgithb/pharmacometrics.git
cd pharmacometrics
bash reproduce.sh
```

`reproduce.sh` restores R packages via `renv`, syncs the Python environment
via `uv`, runs all simulate → fit → diagnose scripts for each project in
order, then calls `quarto render` to build the static site in `docs/`.

Full runtime on 8 cores is approximately 60–90 minutes; the bootstrap
steps (n = 200 replicates, Projects A and B) dominate. Set the environment
variable `BOOTSTRAP_N_BOOT=50` to run a faster smoke-test version.

---

## Repository Structure

```
.
├── R/                        # Shared R helpers (gof_plots, pcvpc, bootstrap_nlmixr, shrinkage_table)
├── python/                   # Shared Python helpers (PBPK ODE builder, plotting utilities)
├── data/simulated/           # Populated by 01_simulate_data.* scripts
├── projects/
│   ├── 01_popPK_nlmixr2/    # Population PK
│   ├── 02_PKPD_emax/        # Sigmoid Emax PK/PD
│   ├── 03_preclinical_to_human/ # FIH allometric scaling
│   ├── 04_PBPK_allometry/   # Minimal PBPK (Python)
│   └── 05_trial_simulation/ # Trial design simulation
├── docs/                     # Rendered HTML (GitHub Pages target)
├── _quarto.yml               # Site configuration
├── reproduce.sh              # End-to-end pipeline script
├── renv.lock                 # Pinned R dependencies
└── pyproject.toml            # Pinned Python dependencies (uv)
```

---

## Dependencies

| Layer | Manager | Version spec |
|-------|---------|--------------|
| R packages | [renv](https://rstudio.github.io/renv/) | See `renv.lock` |
| Python packages | [uv](https://docs.astral.sh/uv/) | See `pyproject.toml` |
| Reports | [Quarto](https://quarto.org) | ≥ 1.5 |

Key R packages: `nlmixr2` (≥ 4.0), `rxode2` (≥ 4.1), `mrgsolve` (≥ 1.6),
`ggplot2` (≥ 4.0), `gt`, `vpc`, `future.apply`.

Key Python packages: `scipy` (≥ 1.13), `pandas` (≥ 2.2), `matplotlib`
(≥ 3.9), `great-tables` (≥ 0.12).

---

## Data

All datasets are generated synthetically by the `01_simulate_data.*` script
in each project directory. Seeds are fixed (`set.seed()` / `np.random.seed()`)
so results are bit-for-bit reproducible. The data-generating (true) parameter
values are documented at the top of each simulation script, enabling
independent verification of parameter recovery.

---

## License

MIT — see [LICENSE](LICENSE).

---

## Author

[Daniela Cruz Moreira]  
[danielacruzmoreira@gmail.com.com]
