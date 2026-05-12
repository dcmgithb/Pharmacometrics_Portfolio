# 02_fit_model.R — Project 03: Multi-species PK fitting and human prediction
# Fits rat and dog separately with nlmixr2 (SAEM), then applies allometric
# scaling to predict human CL and V, and calculates FIH starting dose.

suppressPackageStartupMessages({
  library(nlmixr2)
  library(dplyr)
  library(here)
})

source(here("R", "helpers.R"))

results_dir <- here("projects", "03_preclinical_to_human", "results")
data_dir    <- here("projects", "03_preclinical_to_human", "data")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

rat_data <- read.csv(file.path(data_dir, "rat_pk.csv"))
dog_data <- read.csv(file.path(data_dir, "dog_pk.csv"))

cat("[02_fit_model.R] Data loaded.\n")

# ── IV 1-compartment model (reused for both species) ─────────────────────────
make_iv_model <- function(cl_init, v_init) {
  function() {
    ini({
      log_CL <- log(cl_init)
      log_V  <- log(v_init)
      eta_CL ~ 0.10
      eta_V  ~ 0.05
      prop_err <- 0.20
    })
    model({
      CL <- exp(log_CL) * exp(eta_CL)
      V  <- exp(log_V)  * exp(eta_V)
      d/dt(central) <- -(CL / V) * central
      CONC <- central / V
      CONC ~ prop(prop_err)
    })
  }
}

fit_species <- function(data, cl_init, v_init, label) {
  cat(sprintf("[02_fit_model.R] Fitting %s...\n", label))
  mdl <- make_iv_model(cl_init, v_init)
  fit <- tryCatch(
    nlmixr2(mdl, data, est = "saem",
            control = saemControl(nBurn = 150, nEm = 200, seed = 300, print = 0)),
    error = function(e) {
      warning(label, " fit failed: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(fit))
    cat(sprintf("[02_fit_model.R] %s OFV: %.2f | CL=%.4f | V=%.3f\n",
                label, fit$OBJF,
                exp(fit$fixef["log_CL"]),
                exp(fit$fixef["log_V"])))
  fit
}

# Rat fit: small animal — CL and V both small
fit_rat <- fit_species(rat_data, cl_init = 0.018, v_init = 0.28, label = "Rat")
# Dog fit: large animal — CL and V scaled up
fit_dog <- fit_species(dog_data, cl_init = 1.10,  v_init = 14.0, label = "Dog")

if (!is.null(fit_rat)) saveRDS(fit_rat, file.path(results_dir, "fit_rat.rds"))
if (!is.null(fit_dog)) saveRDS(fit_dog, file.path(results_dir, "fit_dog.rds"))

# ── Allometric scaling ────────────────────────────────────────────────────────
# Fit power law: log(CL) = log(alpha) + beta * log(BW) using the two species.
# With only two data points the exponent cannot be estimated independently —
# standard practice is to fix beta_CL = 0.75 and beta_V = 1.0 (theory-based),
# then estimate the allometric coefficient alpha from the fitted species means.

BW_species <- c(rat = 0.25, dog = 10.0, human = 70.0)

CL_rat_est <- if (!is.null(fit_rat)) exp(fit_rat$fixef["log_CL"]) else 0.020
V_rat_est  <- if (!is.null(fit_rat)) exp(fit_rat$fixef["log_V"])  else 0.30
CL_dog_est <- if (!is.null(fit_dog)) exp(fit_dog$fixef["log_CL"]) else 1.20
V_dog_est  <- if (!is.null(fit_dog)) exp(fit_dog$fixef["log_V"])  else 15.0

# Allometric coefficients: alpha = CL / BW^beta (separately per species, then average)
alpha_CL_rat <- CL_rat_est / BW_species["rat"]^0.75
alpha_CL_dog <- CL_dog_est / BW_species["dog"]^0.75
alpha_V_rat  <- V_rat_est  / BW_species["rat"]^1.00
alpha_V_dog  <- V_dog_est  / BW_species["dog"]^1.00

# Geometric mean of species-specific coefficients
alpha_CL <- exp(mean(log(c(alpha_CL_rat, alpha_CL_dog))))
alpha_V  <- exp(mean(log(c(alpha_V_rat,  alpha_V_dog))))

CL_human <- alpha_CL * BW_species["human"]^0.75
V_human  <- alpha_V  * BW_species["human"]^1.00
t_half   <- log(2) * V_human / CL_human

cat(sprintf(
  "[02_fit_model.R] Human prediction: CL=%.2f L/hr, V=%.1f L, t½=%.1f h\n",
  CL_human, V_human, t_half
))

# ── FIH starting dose (MABEL + NOAEL approach) ────────────────────────────────
# MABEL: minimum anticipated biological effect level
# Assumed in-vitro parameters (stated assumptions, not estimated):
EC50_invitro_nM <- 1.0    # nM (from binding assay, illustrative)
fu_plasma       <- 0.10   # unbound fraction in plasma (protein binding)
MW_kDa          <- 0.400  # molecular weight, kg/mol (400 Da small molecule)

# Convert EC50 to mg/L:  EC50_mgL = EC50_nM * MW (g/mol) / 1e6 (nM→mol/L → mg/L)
MW_gpmol     <- MW_kDa * 1000
EC50_mgL     <- EC50_invitro_nM * MW_gpmol / 1e6  # mg/L

# Plasma concentration at 10% receptor occupancy (MABEL):
# RO = C / (EC50 + C) → C_MABEL = EC50 / (1/0.10 - 1) = EC50 * 0.10 / 0.90
C_mabel_free  <- EC50_mgL * 0.10 / 0.90       # free plasma concentration (mg/L)
C_mabel_total <- C_mabel_free / fu_plasma       # total plasma (mg/L)

# Dose to achieve C_mabel_total at Cmax ≈ C_mabel_total:
# Simple IV: Cmax = Dose / V_human → Dose_MABEL = C_mabel * V_human
Dose_mabel_mg <- C_mabel_total * V_human

# NOAEL-derived HED:
# Assumed: rat 28-day NOAEL = 100 mg/kg; dog NOAEL = 30 mg/kg
NOAEL_rat_mgkg <- 100
NOAEL_dog_mgkg  <- 30

# HED (human equivalent dose): scale by BW^0.75 ratio (body surface area)
HED_rat  <- NOAEL_rat_mgkg * 0.25 * (70 / 0.25)^(1 - 0.75)
HED_dog  <- NOAEL_dog_mgkg * 10.0 * (70 / 10.0)^(1 - 0.75)
HED_most_sensitive <- min(HED_rat, HED_dog)

# Apply 10-fold safety factor (ICH M3(R2) for pharmaceuticals)
Dose_noael_mg <- HED_most_sensitive / 10

# FIH starting dose = lower of MABEL-based and NOAEL-based
Dose_FIH_mg   <- min(Dose_mabel_mg, Dose_noael_mg)
Safety_margin_to_noael <- HED_most_sensitive / Dose_FIH_mg

fih_table <- data.frame(
  Approach         = c("MABEL-based (10% RO)", "NOAEL-based (÷10 safety factor)",
                       "Selected FIH dose"),
  Dose_mg          = round(c(Dose_mabel_mg, Dose_noael_mg, Dose_FIH_mg), 2),
  Safety_margin    = round(c(HED_most_sensitive / Dose_mabel_mg,
                              10,
                              Safety_margin_to_noael), 1),
  Basis            = c(
    sprintf("EC50=%.1f nM, fu=%.2f, V=%.1f L", EC50_invitro_nM, fu_plasma, V_human),
    sprintf("Rat NOAEL=%g mg/kg, HED=%.1f mg, SF=10×", NOAEL_rat_mgkg, HED_most_sensitive),
    "Lower of above"
  )
)
write.csv(fih_table, file.path(results_dir, "fih_dose_table.csv"), row.names = FALSE)

# Human PK prediction table
human_pred <- data.frame(
  Parameter   = c("CL (L/hr)", "V (L)", "t½ (h)", "Alpha_CL", "Alpha_V"),
  Estimate    = round(c(CL_human, V_human, t_half, alpha_CL, alpha_V), 3),
  CI_Lower    = NA_real_,  # populated by 03_diagnostics.R (bootstrap)
  CI_Upper    = NA_real_,
  Notes       = c("Predicted at 70 kg", "Predicted at 70 kg",
                  "Derived: ln2 * V / CL", "Allometric coeff (CL)", "Allometric coeff (V)")
)
write.csv(human_pred, file.path(results_dir, "human_pk_prediction.csv"), row.names = FALSE)

cat(sprintf("[02_fit_model.R] FIH starting dose: %.2f mg\n", Dose_FIH_mg))
cat("[02_fit_model.R] Done. Run 03_diagnostics.R next.\n")
