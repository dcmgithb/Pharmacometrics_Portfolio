# 02_fit_model.R — Project 02: Sequential PK→PD estimation (sigmoid Emax)
# Step 1: Fit 1-cmt PK with nlmixr2 SAEM.
# Step 2: Fix PK parameters; fit sigmoid Emax PD model.
# Step 3: Compare Cmax, AUC, Cavg as exposure metrics (AIC).

suppressPackageStartupMessages({
  library(nlmixr2)
  library(dplyr)
  library(here)
})

source(here("R", "helpers.R"))

results_dir <- here("projects", "02_PKPD_emax", "results")
data_dir    <- here("projects", "02_PKPD_emax", "data")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

pkpd_data    <- read.csv(file.path(data_dir, "pkpd_data.csv"))
exposure_data <- read.csv(file.path(data_dir, "exposure_metrics.csv"))

cat("[02_fit_model.R] Data loaded:", nrow(pkpd_data), "rows\n")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1 — PK sub-model (1-compartment oral)
# ═══════════════════════════════════════════════════════════════════════════════
pk_only <- pkpd_data[pkpd_data$CMT == 1, ]

pk_model <- function() {
  ini({
    log_CL <- log(2.5)   # true 3.0, init 2.5
    log_V  <- log(18)    # true 20, init 18
    log_Ka <- log(0.7)   # true 0.8, init 0.7
    eta_CL ~ 0.10
    eta_V  ~ 0.05
    prop_err <- 0.20
  })
  model({
    CL <- exp(log_CL) * exp(eta_CL)
    V  <- exp(log_V)  * exp(eta_V)
    Ka <- exp(log_Ka)
    d/dt(depot)   <- -Ka * depot
    d/dt(central) <- Ka * depot - (CL/V) * central
    CONC <- central / V
    CONC ~ prop(prop_err)
  })
}

cat("[02_fit_model.R] Fitting PK sub-model (SAEM)...\n")
fit_pk <- nlmixr2(pk_model, pk_only,
                  est = "saem",
                  control = saemControl(nBurn = 200, nEm = 300, seed = 2001, print = 50))
saveRDS(fit_pk, file.path(results_dir, "fit_pk.rds"))

cat("[02_fit_model.R] PK OFV:", fit_pk$OBJF, "\n")

# Extract individual AUC via IPRED from the PK fit (model-based NCA)
pk_fit_df <- as.data.frame(fit_pk)
pk_obs    <- pk_fit_df[pk_fit_df$EVID == 0, ]

# Compute individual AUC0-24 from IPRED using trapezoidal rule
trapz <- function(x, y) {
  idx <- order(x); sum(diff(x[idx]) * (y[idx[-length(y)]] + y[idx[-1]]) / 2)
}

indiv_metrics <- pk_obs %>%
  group_by(ID) %>%
  summarise(
    AUC_fit = trapz(TIME, IPRED),
    Cmax_fit = max(IPRED),
    Cavg_fit = mean(IPRED[TIME <= 24]),
    .groups = "drop"
  )

# Merge with PD observations
pd_obs <- pkpd_data[pkpd_data$CMT == 2, c("ID","DV","DOSE")]
pd_fit_data <- left_join(pd_obs, indiv_metrics, by = "ID") %>%
  left_join(exposure_data[, c("ID","AUC0_24","Cmax","Cavg")], by = "ID", suffix = c("","_true"))

write.csv(pd_fit_data, file.path(results_dir, "pd_fit_data.csv"), row.names = FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2 — PD sub-model: Emax with three exposure metrics compared
# ═══════════════════════════════════════════════════════════════════════════════
# Use FOCE-I (focei) for PD sub-model: sparser PD data makes SAEM less stable,
# and FOCE-I is more reliable for models where individual predictions are
# analytically differentiable (algebraic Emax, no ODEs).

fit_emax_with_metric <- function(metric_col, data) {
  data$EXPOSURE <- data[[metric_col]]

  pd_model <- function() {
    ini({
      log_E0   <- log(95)   # true 100
      log_Emax <- log(75)   # true 80
      log_EC50 <- log(2.0)  # true 2.5
      log_gamma<- log(1.5)  # true 1.8
      eta_E0   ~ 0.05
      eta_Emax ~ 0.10
      eta_EC50 ~ 0.20
      add_err  <- 10
    })
    model({
      E0    <- exp(log_E0)   * exp(eta_E0)
      Emax  <- exp(log_Emax) * exp(eta_Emax)
      EC50  <- exp(log_EC50) * exp(eta_EC50)
      gamma <- exp(log_gamma)
      EFFECT <- E0 + (Emax * EXPOSURE^gamma) / (EC50^gamma + EXPOSURE^gamma)
      EFFECT ~ add(add_err)
    })
  }

  tryCatch(
    nlmixr2(pd_model, data,
            est = "focei",
            control = foceiControl(maxEval = 10000, print = 0)),
    error = function(e) {
      warning(metric_col, ": ", conditionMessage(e))
      NULL
    }
  )
}

cat("[02_fit_model.R] Fitting PD model with AUC metric...\n")
fit_auc  <- fit_emax_with_metric("AUC_fit",  pd_fit_data)
cat("[02_fit_model.R] Fitting PD model with Cmax metric...\n")
fit_cmax <- fit_emax_with_metric("Cmax_fit", pd_fit_data)
cat("[02_fit_model.R] Fitting PD model with Cavg metric...\n")
fit_cavg <- fit_emax_with_metric("Cavg_fit", pd_fit_data)

# ── Metric comparison table ────────────────────────────────────────────────────
fits_list <- list(AUC = fit_auc, Cmax = fit_cmax, Cavg = fit_cavg)
metric_df <- lapply(names(fits_list), function(nm) {
  f <- fits_list[[nm]]
  if (is.null(f)) return(data.frame(Metric=nm, OFV=NA, AIC=NA, BIC=NA))
  data.frame(Metric = nm,
             OFV    = round(f$OBJF, 2),
             AIC    = round(f$AIC,  2),
             BIC    = round(f$BIC,  2))
}) %>% do.call(rbind)

metric_df <- metric_df[order(metric_df$AIC), ]
write.csv(metric_df, file.path(results_dir, "metric_comparison.csv"), row.names = FALSE)
cat("[02_fit_model.R] Metric comparison:\n"); print(metric_df)

# Keep the best fit (lowest AIC)
best_metric <- metric_df$Metric[which.min(metric_df$AIC)]
fit_pd_best <- fits_list[[best_metric]]
cat("[02_fit_model.R] Best metric:", best_metric, "\n")

if (!is.null(fit_pd_best)) {
  saveRDS(fit_pd_best, file.path(results_dir, "fit_pd_best.rds"))

  # PD parameter table
  pd_pars <- data.frame(
    Parameter = names(fit_pd_best$fixef),
    Estimate  = round(unname(fit_pd_best$fixef), 4),
    AIC       = round(fit_pd_best$AIC, 2),
    Metric    = best_metric
  )
  write.csv(pd_pars, file.path(results_dir, "pd_parameters.csv"), row.names = FALSE)
}

# ── Phase 2 dose recommendation ───────────────────────────────────────────────
# Find dose range achieving 50–80% of Emax at population median
if (!is.null(fit_pd_best)) {
  theta_pd <- fit_pd_best$fixef

  E0_est   <- exp(theta_pd["log_E0"])
  Emax_est <- exp(theta_pd["log_Emax"])
  EC50_est <- exp(theta_pd["log_EC50"])
  gamma_est<- exp(theta_pd["log_gamma"])

  # Dose-to-AUC conversion from PK fit
  CL_est   <- exp(fit_pk$fixef["log_CL"])
  F_oral   <- 0.85  # assumed bioavailability (not estimated)
  dose_seq <- seq(10, 600, by = 5)
  AUC_seq  <- (dose_seq * F_oral) / CL_est   # steady-state AUC ≈ dose*F/CL

  EFFECT_seq <- E0_est + (Emax_est * AUC_seq^gamma_est) /
                          (EC50_est^gamma_est + AUC_seq^gamma_est)
  pct_emax   <- (EFFECT_seq - E0_est) / Emax_est

  phase2_df <- data.frame(
    Dose_mg     = dose_seq,
    AUC_pred    = round(AUC_seq, 2),
    Effect_pred = round(EFFECT_seq, 1),
    PctEmax     = round(pct_emax * 100, 1)
  )

  rec_df <- phase2_df[phase2_df$PctEmax >= 50 & phase2_df$PctEmax <= 80, ]
  rec_summary <- data.frame(
    Dose_min_mg   = min(rec_df$Dose_mg),
    Dose_max_mg   = max(rec_df$Dose_mg),
    PctEmax_min   = min(rec_df$PctEmax),
    PctEmax_max   = max(rec_df$PctEmax),
    Metric_used   = best_metric,
    EC50_estimate = round(EC50_est, 3),
    EC50_true     = 2.5
  )
  write.csv(rec_summary, file.path(results_dir, "phase2_recommendation.csv"), row.names = FALSE)
  write.csv(phase2_df,   file.path(results_dir, "er_curve.csv"), row.names = FALSE)
  cat("[02_fit_model.R] Phase 2 recommendation:", rec_summary$Dose_min_mg,
      "–", rec_summary$Dose_max_mg, "mg\n")
}

cat("[02_fit_model.R] Done. Run 03_diagnostics.R next.\n")
