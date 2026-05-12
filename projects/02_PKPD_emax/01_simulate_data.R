# 01_simulate_data.R — Project 02: PK/PD simulation (1-cmt PK + sigmoid Emax)
#
# ═══════════════════════════════════════════════════════════════════════════════
# TRUE PARAMETER VALUES
# ═══════════════════════════════════════════════════════════════════════════════
# PK (1-compartment oral):
#   CL = 3.00 L/hr,  V = 20.0 L,  Ka = 0.80 1/hr
#   IIV: omega_CL = 0.09 (CV≈30%), omega_V = 0.04 (CV≈20%)
#   RUV: prop_pk = 0.15 (CV≈15%)
#
# PD (direct sigmoid Emax, driven by plasma AUC metric):
#   E0    = 100.0 (baseline; e.g. pain score, 0-200 VAS)
#   Emax  =  80.0 (maximum achievable effect)
#   EC50  =   2.5 mg·h/L  (AUC0-24 at half-maximum effect)
#   gamma =   1.8 (Hill coefficient — moderate sigmoidicity)
#   IIV: omega_E0 = 0.04 (CV≈20%), omega_Emax = 0.09 (CV≈30%)
#         omega_EC50 = 0.16 (CV≈40%) — EC50 is the most variable PD parameter
#   RUV: sigma_pd_add = 8.0 (absolute units)
#
# Design:
#   4 dose groups: 50, 100, 200, 400 mg oral QD; n=20 per group (N=80)
#   PK sampling: 0.5, 1, 2, 4, 8, 24 h post-dose (Day 1 single-dose)
#   PD: one biomarker observation per subject at 24 h (trough effect)
#   AUC0-24 computed from individual IPRED (NCA trapezoidal)
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(rxode2)
  library(dplyr)
})

set.seed(2025)

# ── True parameters ───────────────────────────────────────────────────────────
CL_true  <- 3.00;  V_true   <- 20.0;  Ka_true  <- 0.80
om_CL    <- 0.09;  om_V     <- 0.04
sig_pk   <- 0.15   # proportional CV

E0_true   <- 100.0; Emax_true <- 80.0; EC50_true <- 2.5; gamma_true <- 1.8
om_E0     <- 0.04;  om_Emax  <- 0.09; om_EC50   <- 0.16
sig_pd    <- 8.0   # additive (units of the biomarker scale)

DOSES     <- c(50, 100, 200, 400)  # mg
N_PER_ARM <- 20
PK_TIMES  <- c(0.5, 1, 2, 4, 8, 24)  # hours post-dose

# ── rxode2 PK model ───────────────────────────────────────────────────────────
pk_model <- rxode2({
  d/dt(depot)   <- -Ka * depot
  d/dt(central) <- Ka * depot - (CL/V) * central
  CONC          <- central / V
})

# ── Individual parameters ─────────────────────────────────────────────────────
N_total <- length(DOSES) * N_PER_ARM

DOSE_vec <- rep(DOSES, each = N_PER_ARM)

eta_CL  <- rnorm(N_total, 0, sqrt(om_CL))
eta_V   <- rnorm(N_total, 0, sqrt(om_V))
eta_E0  <- rnorm(N_total, 0, sqrt(om_E0))
eta_Em  <- rnorm(N_total, 0, sqrt(om_Emax))
eta_EC  <- rnorm(N_total, 0, sqrt(om_EC50))

CL_i   <- CL_true * exp(eta_CL)
V_i    <- V_true  * exp(eta_V)
Ka_i   <- Ka_true  # no IIV on Ka for simplicity
E0_i   <- E0_true  * exp(eta_E0)
Emax_i <- Emax_true * exp(eta_Em)
EC50_i <- EC50_true * exp(eta_EC)

# ── Simulate PK ───────────────────────────────────────────────────────────────
pk_records  <- list()
pd_records  <- list()

for (i in seq_len(N_total)) {
  dose <- DOSE_vec[i]
  ev   <- et(amt = dose, time = 0, cmt = "depot") %>% et(PK_TIMES)

  sol <- rxSolve(pk_model,
                 c(Ka = Ka_i, CL = CL_i[i], V = V_i[i]),
                 ev,
                 returnType = "data.frame",
                 seed = 2025 + i)

  pk_obs <- sol[sol$time %in% PK_TIMES, ]
  pk_obs$DV  <- pk_obs$CONC * (1 + rnorm(nrow(pk_obs), 0, sig_pk))
  pk_obs$DV  <- pmax(pk_obs$DV, 0)
  pk_obs$ID  <- i
  pk_obs$CMT <- 1
  pk_obs$EVID<- 0
  pk_obs$MDV <- 0
  pk_obs$AMT <- 0
  pk_obs$DOSE<- dose
  pk_records[[i]] <- pk_obs[, c("ID","time","DV","AMT","EVID","MDV","CMT","CONC","DOSE")]

  # Individual AUC0-24 from trapezoidal rule on individual predicted concentrations
  AUC_i <- trapz(sol$time, sol$CONC)  # mg·h/L

  # PD response driven by AUC (true exposure metric in this simulation)
  EFFECT_i <- E0_i[i] + (Emax_i[i] * AUC_i^gamma_true) /
                         (EC50_i[i]^gamma_true + AUC_i^gamma_true)
  DV_pd    <- EFFECT_i + rnorm(1, 0, sig_pd)

  # Also compute Cmax and Cavg for the metric comparison analysis
  Cmax_i <- max(sol$CONC)
  Cavg_i <- mean(sol$CONC[sol$time <= 24])

  pd_records[[i]] <- data.frame(
    ID   = i, TIME = 24, DV = DV_pd, AMT = 0, EVID = 0, MDV = 0,
    CMT  = 2, DOSE = dose,
    AUC0_24 = AUC_i, Cmax = Cmax_i, Cavg = Cavg_i,
    EFFECT_true = EFFECT_i
  )
}

# ── Assemble NONMEM-style dataset ─────────────────────────────────────────────
pk_df <- do.call(rbind, pk_records)
names(pk_df)[names(pk_df) == "time"] <- "TIME"

pd_df <- do.call(rbind, pd_records)

# Merge PK and PD with dose records
dose_rows <- data.frame(
  ID = seq_len(N_total), TIME = 0, DV = NA, AMT = DOSE_vec,
  EVID = 1, MDV = 1, CMT = 1, DOSE = DOSE_vec,
  CONC = NA, AUC0_24 = NA, Cmax = NA, Cavg = NA, EFFECT_true = NA
)

pkpd_data <- dplyr::bind_rows(
  dplyr::select(pk_df, ID, TIME, DV, AMT, EVID, MDV, CMT, DOSE),
  dplyr::select(pd_df, ID, TIME, DV, AMT, EVID, MDV, CMT, DOSE)
) %>% arrange(ID, TIME)

# Save exposure metrics separately for the metric comparison analysis
exposure_data <- pd_df[, c("ID", "DOSE", "AUC0_24", "Cmax", "Cavg",
                            "DV", "EFFECT_true")]
names(exposure_data)[names(exposure_data) == "DV"] <- "DV_PD"

out_dir <- "data"
dir.create(out_dir, showWarnings = FALSE)
write.csv(pkpd_data,     file.path(out_dir, "pkpd_data.csv"),     row.names = FALSE)
write.csv(exposure_data, file.path(out_dir, "exposure_metrics.csv"), row.names = FALSE)

cat(sprintf(
  "[01_simulate_data.R] Done.\n  N subjects: %d | Dose groups: %s\n  PK obs: %d | PD obs: %d\n",
  N_total, paste(DOSES, collapse="/"), nrow(pk_df), nrow(pd_df)
))

# Helper: trapezoidal integration (not loaded by default in base R)
trapz <- function(x, y) {
  idx <- order(x)
  sum(diff(x[idx]) * (y[idx][-length(y)] + y[idx][-1]) / 2)
}
