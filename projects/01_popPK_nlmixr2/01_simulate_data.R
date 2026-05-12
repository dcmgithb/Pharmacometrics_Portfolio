# 01_simulate_data.R — Project 01: Population PK simulation
# Generates a synthetic sparse-sampling PK dataset from a known 2-compartment
# oral model with allometric weight on CL and V.
#
# ═══════════════════════════════════════════════════════════════════════════════
# TRUE PARAMETER VALUES (data-generating model)
# ═══════════════════════════════════════════════════════════════════════════════
# CL_ref  = 5.00 L/hr  (clearance at 70 kg reference body weight)
# V2_ref  = 50.0 L     (central volume at reference BW)
# V3_ref  = 100. L     (peripheral volume; not weight-scaled — parsimony)
# Q_ref   = 3.00 L/hr  (intercompartmental clearance)
# Ka      = 1.20 1/hr  (first-order absorption rate constant)
# theta_CL_BW = 0.75   (allometric exponent on CL)
# theta_V_BW  = 1.00   (allometric exponent on V2)
# omega_CL    = 0.09   (variance; CV ≈ 30%)
# omega_V2    = 0.0625 (variance; CV ≈ 25%)
# omega_Ka    = 0.16   (variance; CV ≈ 40%)
# sigma_prop  = 0.0400 (proportional RUV variance; CV ≈ 20%)
# sigma_add   = 0.0100 (additive RUV variance; SD ≈ 0.1 ng/mL)
# ═══════════════════════════════════════════════════════════════════════════════
# N subjects = 120, body weight ~ N(75, 15²) truncated [40, 120] kg
# Dose: 100 mg oral QD for 14 days (steady-state trough available)
# Sampling: 3–5 sparse time points per subject (occasions: Day 1 + Day 14)
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(rxode2)
  library(dplyr)
  library(tidyr)
})

set.seed(2024)  # Fixed seed: ensures bit-for-bit reproducible output

# ── True parameters ───────────────────────────────────────────────────────────
CL_ref       <- 5.00
V2_ref       <- 50.0
V3_ref       <- 100.0
Q_ref        <- 3.00
Ka_true      <- 1.20
theta_CL_BW  <- 0.75
theta_V_BW   <- 1.00
omega_CL     <- 0.09
omega_V2     <- 0.0625
omega_Ka     <- 0.16
sigma_prop   <- 0.04
sigma_add    <- 0.01

N_subj       <- 120
DOSE_MG      <- 100   # mg

# ── Covariate generation ──────────────────────────────────────────────────────
# Truncated normal for body weight: realistic clinical range
rtnorm <- function(n, mean, sd, lo, hi) {
  x <- rnorm(n, mean, sd)
  x[x < lo] <- lo + abs(rnorm(sum(x < lo), 0, sd / 3))
  x[x > hi] <- hi - abs(rnorm(sum(x > hi), 0, sd / 3))
  pmin(pmax(x, lo), hi)
}

BW  <- rtnorm(N_subj, mean = 75, sd = 15, lo = 40, hi = 120)
AGE <- round(runif(N_subj, 18, 65))
SEX <- rbinom(N_subj, 1, 0.50)   # 1 = female

# ── Individual PK parameters (lognormal IIV) ─────────────────────────────────
eta_CL <- rnorm(N_subj, 0, sqrt(omega_CL))
eta_V2 <- rnorm(N_subj, 0, sqrt(omega_V2))
eta_Ka <- rnorm(N_subj, 0, sqrt(omega_Ka))

CL_i <- CL_ref * (BW / 70)^theta_CL_BW * exp(eta_CL)
V2_i <- V2_ref * (BW / 70)^theta_V_BW  * exp(eta_V2)
V3_i <- rep(V3_ref, N_subj)   # no covariate on peripheral volume
Q_i  <- rep(Q_ref,  N_subj)
Ka_i <- Ka_true * exp(eta_Ka)

# ── Sampling design ───────────────────────────────────────────────────────────
# Sparse design: Day 1 (0, 1, 4 h post-dose) + Day 14 (pre-dose, 2, 24 h)
# Each subject gets 3–5 of these time points (random dropout)
samp_day1 <- c(0.5, 1, 2, 4, 8)     # post Day-1 dose
samp_day14<- c(0, 1, 3, 6, 24) + 13*24  # Day 14 (hours from first dose)

# ── rxode2 model definition ───────────────────────────────────────────────────
pk2cmt_model <- rxode2({
  d/dt(depot)      <- -Ka * depot
  d/dt(central)    <- Ka * depot - (CL/V2) * central - (Q/V2) * central + (Q/V3) * periph
  d/dt(periph)     <- (Q/V2) * central - (Q/V3) * periph
  CONC             <- central / V2
  DV               <- CONC * (1 + eps_prop) + eps_add
})

# ── Simulate each subject ─────────────────────────────────────────────────────
all_data <- lapply(seq_len(N_subj), function(i) {
  # Sparse sample selection (3-5 obs per occasion)
  n_d1  <- sample(2:4, 1)
  n_d14 <- sample(2:4, 1)
  obs_times <- c(
    sort(sample(samp_day1,  n_d1,  replace = FALSE)),
    sort(sample(samp_day14, n_d14, replace = FALSE))
  )

  # Event table: doses every 24h for 14 days + observation records
  ev <- et(amt  = DOSE_MG,
            time = seq(0, 13*24, by = 24),
            cmt  = "depot",
            addl = 0) %>%
        et(obs_times)

  # Simulate with individual parameters + RUV
  theta_i <- c(Ka = Ka_i[i], CL = CL_i[i], V2 = V2_i[i],
               V3 = V3_i[i], Q  = Q_i[i])

  sim_i <- rxSolve(pk2cmt_model, theta_i, ev,
                   omega  = NULL,  # IIV already incorporated above
                   sigma  = matrix(c(sigma_prop, 0, 0, sigma_add), 2, 2),
                   sigmaLower = 0,
                   returnType = "data.frame",
                   seed    = 2024 + i)

  # Keep only observation rows and format as NONMEM-style dataset
  obs_rows <- sim_i[sim_i$time %in% obs_times, ]
  obs_rows$DV[obs_rows$DV < 0] <- 0  # concentrations cannot be negative

  dose_rows <- data.frame(
    ID   = i, TIME = seq(0, 13*24, by = 24),
    DV   = NA_real_, AMT = DOSE_MG, EVID = 1, MDV = 1,
    CMT  = 1, BW = BW[i], AGE = AGE[i], SEX = SEX[i], DOSE = DOSE_MG
  )
  obs_out <- data.frame(
    ID   = i, TIME = obs_rows$time,
    DV   = pmax(obs_rows$DV, 0), AMT = 0, EVID = 0, MDV = 0,
    CMT  = 2, BW = BW[i], AGE = AGE[i], SEX = SEX[i], DOSE = DOSE_MG
  )
  rbind(dose_rows, obs_out)
})

pk_data <- do.call(rbind, all_data)
pk_data  <- pk_data[order(pk_data$ID, pk_data$TIME), ]

# Weight tertile covariate (used for VPC stratification)
pk_data$BW_tertile <- cut(pk_data$BW,
                           breaks = quantile(pk_data$BW, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                           labels = c("Low (<60kg)", "Mid (60-85kg)", "High (>85kg)"),
                           include.lowest = TRUE)

# ── Rich-sampling validation set (20 subjects, intensive PK) ─────────────────
set.seed(20240)
N_rich    <- 20
BW_rich   <- rtnorm(N_rich, 70, 12, 45, 110)
eta_CL_r  <- rnorm(N_rich, 0, sqrt(omega_CL))
eta_V2_r  <- rnorm(N_rich, 0, sqrt(omega_V2))
eta_Ka_r  <- rnorm(N_rich, 0, sqrt(omega_Ka))
CL_r <- CL_ref * (BW_rich/70)^theta_CL_BW * exp(eta_CL_r)
V2_r <- V2_ref * (BW_rich/70)^theta_V_BW  * exp(eta_V2_r)
Ka_r <- Ka_true * exp(eta_Ka_r)

rich_times <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 24)  # single-dose Day 1

rich_data <- lapply(seq_len(N_rich), function(i) {
  ev <- et(amt = DOSE_MG, time = 0, cmt = "depot") %>% et(rich_times)
  theta_i <- c(Ka = Ka_r[i], CL = CL_r[i], V2 = V2_ref,
               V3 = V3_ref,  Q  = Q_ref)
  sim_i   <- rxSolve(pk2cmt_model, theta_i, ev,
                     sigma = matrix(c(sigma_prop, 0, 0, sigma_add), 2, 2),
                     returnType = "data.frame", seed = 9000 + i)
  obs_rows <- sim_i[sim_i$time %in% rich_times, ]
  data.frame(ID = N_subj + i, TIME = obs_rows$time,
             DV = pmax(obs_rows$DV, 0), AMT = 0, EVID = 0, MDV = 0,
             CMT = 2, BW = BW_rich[i], AGE = 40, SEX = 0,
             DOSE = DOSE_MG, BW_tertile = "Rich")
})
pk_data_rich <- do.call(rbind, rich_data)

# ── Save datasets ─────────────────────────────────────────────────────────────
out_dir <- file.path("data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(pk_data,      file.path(out_dir, "pk_data.csv"),      row.names = FALSE)
write.csv(pk_data_rich, file.path(out_dir, "pk_data_rich.csv"), row.names = FALSE)

cat(sprintf(
  "[01_simulate_data.R] Done.\n  Sparse dataset:  %d subjects, %d observations\n  Rich dataset:    %d subjects, %d observations\n  Saved to data/\n",
  length(unique(pk_data$ID[pk_data$EVID==0])),
  sum(pk_data$EVID == 0),
  N_rich,
  nrow(pk_data_rich)
))
