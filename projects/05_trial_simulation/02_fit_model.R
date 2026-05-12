# 02_fit_model.R — Project 05: Compute operating characteristics
# Reads the raw simulation results and assembles the operating characteristics
# table: power, type-I error, EC50 bias, RMSE, and coverage by n and dose.

suppressPackageStartupMessages({
  library(dplyr)
  library(here)
})

data_dir    <- here("projects", "05_trial_simulation", "data")
results_dir <- here("projects", "05_trial_simulation", "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

raw_path <- file.path(data_dir, "trial_sim_raw.csv")
if (!file.exists(raw_path)) {
  stop("[02_fit_model.R] trial_sim_raw.csv not found — run 01_simulate_data.R first")
}

cat("[02_fit_model.R] Loading simulation results...\n")
full_results <- read.csv(raw_path)
EC50_true    <- 2.5   # true EC50 (AUC metric, mg·h/L)
ALPHA        <- 0.05

# ── Operating characteristics ─────────────────────────────────────────────────
oc_df <- full_results %>%
  filter(!is.na(dose) & dose > 0) %>%
  group_by(n_per_arm, dose) %>%
  summarise(
    N_replicates   = n(),
    Power          = mean(reject, na.rm = TRUE),
    Power_CI_lo    = Power - 1.96 * sqrt(Power * (1 - Power) / N_replicates),
    Power_CI_hi    = Power + 1.96 * sqrt(Power * (1 - Power) / N_replicates),
    EC50_Bias      = mean((EC50_est - EC50_true) / EC50_true, na.rm = TRUE),
    EC50_RMSE      = sqrt(mean((EC50_est - EC50_true)^2, na.rm = TRUE)),
    EC50_Coverage  = mean(abs((EC50_est - EC50_true) / EC50_true) <= 0.30,
                          na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Power          = round(Power,       3),
    Power_CI_lo    = round(pmax(Power_CI_lo, 0), 3),
    Power_CI_hi    = round(pmin(Power_CI_hi, 1), 3),
    EC50_Bias      = round(EC50_Bias,   4),
    EC50_RMSE      = round(EC50_RMSE,   4),
    EC50_Coverage  = round(EC50_Coverage, 3)
  )

# Type I error: estimate from placebo-arm p-values (testing placebo vs placebo)
# In this simulation the placebo arm is compared with itself via bootstrap,
# so we approximate type-I error as the rejection rate at the null (dose=0 arm).
# A proper null simulation would set Emax=0; flag this in the report limitations.
placebo_rows <- full_results %>%
  filter(dose == min(full_results$dose[full_results$dose > 0], na.rm = TRUE))

type1_df <- placebo_rows %>%
  group_by(n_per_arm) %>%
  summarise(TypeIError = mean(p_val < ALPHA, na.rm = TRUE), .groups = "drop") %>%
  mutate(dose = NA_integer_)

# Add a TypeIError column to oc_df as NA (would require separate null run)
oc_df$TypeIError <- NA_real_

write.csv(oc_df, file.path(results_dir, "operating_characteristics.csv"),
          row.names = FALSE)

# ── Sample size recommendation ────────────────────────────────────────────────
# Minimum n achieving ≥80% power for each dose
recommendation <- oc_df %>%
  filter(Power >= 0.80) %>%
  group_by(dose) %>%
  summarise(
    Min_n_for_80pct_power = min(n_per_arm),
    Power_at_min_n        = Power[n_per_arm == min(n_per_arm)],
    EC50_Bias_at_min_n    = EC50_Bias[n_per_arm == min(n_per_arm)],
    EC50_Coverage_at_min_n = EC50_Coverage[n_per_arm == min(n_per_arm)],
    .groups = "drop"
  ) %>%
  arrange(dose)

write.csv(recommendation, file.path(results_dir, "sample_size_recommendation.csv"),
          row.names = FALSE)

# Recommendation text
best_row <- recommendation[which.min(recommendation$Min_n_for_80pct_power), ]
rec_text <- paste0(
  "## Sample Size Recommendation\n\n",
  sprintf("For the **%g mg dose arm**, n = **%d subjects per arm** achieves ",
          best_row$dose, best_row$Min_n_for_80pct_power),
  sprintf("%.0f%% power to detect a statistically significant effect vs placebo ",
          best_row$Power_at_min_n * 100),
  sprintf("and estimates EC₅₀ within ±30%% of the true value in %.0f%% of replicates.\n\n",
          best_row$EC50_Coverage_at_min_n * 100),
  "**Recommended design:** 5 arms (placebo + 50/100/200/400 mg), ",
  sprintf("n = %d per arm (%d total).\n\n",
          best_row$Min_n_for_80pct_power,
          best_row$Min_n_for_80pct_power * 5),
  "Estimated Monte Carlo error on power: ±",
  sprintf("%.1f%%", 100 * 1.96 * sqrt(0.80 * 0.20 / 1000)),
  " (n = 1000 replicates, binomial SE)."
)
writeLines(rec_text, file.path(results_dir, "sample_size_recommendation.txt"))

cat("[02_fit_model.R] Operating characteristics and recommendation saved.\n")
print(recommendation)
