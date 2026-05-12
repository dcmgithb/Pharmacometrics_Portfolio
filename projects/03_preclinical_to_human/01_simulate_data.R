# 01_simulate_data.R — Project 03: Preclinical PK simulation (rat + dog)
#
# ═══════════════════════════════════════════════════════════════════════════════
# TRUE PARAMETER VALUES
# ═══════════════════════════════════════════════════════════════════════════════
# Rat (BW = 0.25 kg):
#   CL_rat = 0.020 L/hr,  V_rat = 0.30 L
#   (allometric: CL = alpha_CL * BW^0.75 → alpha_CL = 0.020 / 0.25^0.75 = 0.0669)
#
# Dog (BW = 10.0 kg):
#   CL_dog = 1.20 L/hr,   V_dog = 15.0 L
#   (allometric: CL = alpha_CL * 10^0.75 = 0.0669 * 5.62 = 0.376... ≠ 1.20 intentionally
#    — the dog value is set to create slight inter-species scatter, realistic noise)
#
# Human prediction (from allometric fit, not assumed):
#   CL_human_pred ≈ alpha_CL_fit * 70^0.75
#   V_human_pred  ≈ alpha_V_fit  * 70^1.00
#
# IIV within species: omega_CL = 0.09 (CV≈30%), omega_V = 0.04 (CV≈20%)
# RUV: proportional, sigma_prop = 0.15
#
# Design:
#   Rat: n=6, IV bolus 1 mg/kg, rich sampling 0–24 h (9 timepoints)
#   Dog: n=4, IV bolus 1 mg/kg, rich sampling 0–48 h (11 timepoints)
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(rxode2)
  library(dplyr)
})

set.seed(2026)

# ── True parameters ───────────────────────────────────────────────────────────
# Species physiological parameters
species_params <- list(
  rat = list(BW = 0.25, CL = 0.020, V = 0.30,
             dose_mgkg = 1, times = c(0.083, 0.25, 0.5, 1, 2, 4, 8, 12, 24)),
  dog = list(BW = 10.0,  CL = 1.20,  V = 15.0,
             dose_mgkg = 1, times = c(0.25, 0.5, 1, 2, 4, 8, 12, 24, 36, 48))
)
n_animals <- list(rat = 6, dog = 4)
omega_CL   <- 0.09
omega_V    <- 0.04
sigma_prop <- 0.15

# ── rxode2 1-cmt IV model ─────────────────────────────────────────────────────
iv_model <- rxode2({
  d/dt(central) <- -(CL / V) * central
  CONC          <- central / V
})

# ── Simulate per-species ──────────────────────────────────────────────────────
simulate_species <- function(sp_name, n, params, seed_offset = 0) {
  bw     <- params$BW
  cl_ref <- params$CL
  v_ref  <- params$V
  dose   <- params$dose_mgkg * bw  # mg (mg/kg × BW)
  times  <- params$times

  lapply(seq_len(n), function(i) {
    set.seed(2026 + seed_offset + i)
    CL_i <- cl_ref * exp(rnorm(1, 0, sqrt(omega_CL)))
    V_i  <- v_ref  * exp(rnorm(1, 0, sqrt(omega_V)))

    ev  <- et(amt = dose, time = 0, cmt = "central") %>% et(times)
    sol <- rxSolve(iv_model, c(CL = CL_i, V = V_i), ev,
                   returnType = "data.frame", seed = 100 * i)

    obs <- sol[sol$time %in% times, ]
    obs$DV     <- pmax(obs$CONC * (1 + rnorm(nrow(obs), 0, sigma_prop)), 0)
    obs$ID     <- i
    obs$SPECIES <- sp_name
    obs$BW     <- bw
    obs$DOSE   <- dose
    obs$CL_true <- CL_i
    obs$V_true  <- V_i
    obs$EVID   <- 0
    obs$MDV    <- 0
    obs$CMT    <- 1
    obs$AMT    <- 0

    obs[, c("ID","SPECIES","BW","time","DV","AMT","EVID","MDV","CMT",
            "DOSE","CL_true","V_true")]
  }) %>% do.call(rbind)
}

rat_data <- simulate_species("rat", n_animals$rat, species_params$rat, seed_offset = 0)
dog_data <- simulate_species("dog", n_animals$dog, species_params$dog, seed_offset = 100)

# Rename time column for consistency
names(rat_data)[names(rat_data) == "time"] <- "TIME"
names(dog_data)[names(dog_data) == "time"] <- "TIME"

# Add dose records
make_dose_rows <- function(data, sp_name, params) {
  ids <- unique(data$ID)
  lapply(ids, function(id) {
    bw   <- params$BW
    dose <- params$dose_mgkg * bw
    data.frame(ID = id, SPECIES = sp_name, BW = bw, TIME = 0,
               DV = NA, AMT = dose, EVID = 1, MDV = 1, CMT = 1,
               DOSE = dose, CL_true = NA, V_true = NA)
  }) %>% do.call(rbind)
}

rat_doses <- make_dose_rows(rat_data, "rat", species_params$rat)
dog_doses <- make_dose_rows(dog_data, "dog", species_params$dog)

rat_full <- rbind(rat_doses, rat_data) %>% arrange(ID, TIME)
dog_full <- rbind(dog_doses, dog_data) %>% arrange(ID, TIME)

out_dir <- "data"
dir.create(out_dir, showWarnings = FALSE)
write.csv(rat_full, file.path(out_dir, "rat_pk.csv"), row.names = FALSE)
write.csv(dog_full, file.path(out_dir, "dog_pk.csv"), row.names = FALSE)

cat(sprintf(
  "[01_simulate_data.R] Done.\n  Rat: %d animals, %d obs | Dog: %d animals, %d obs\n",
  n_animals$rat, sum(rat_data$EVID == 0),
  n_animals$dog, sum(dog_data$EVID == 0)
))
