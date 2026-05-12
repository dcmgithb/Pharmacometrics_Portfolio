# 03_diagnostics.R — Project 05: Power curve, EC50 precision, OC table plots

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(gt)
  library(here)
})

source(here("R", "helpers.R"))

results_dir <- here("projects", "05_trial_simulation", "results")
data_dir    <- here("projects", "05_trial_simulation", "data")

oc_path <- file.path(results_dir, "operating_characteristics.csv")
if (!file.exists(oc_path)) stop("Run 02_fit_model.R first (operating_characteristics.csv missing)")

oc_df  <- read.csv(oc_path)
raw_df <- read.csv(file.path(data_dir, "trial_sim_raw.csv"))

cat("[03_diagnostics.R] Generating trial simulation diagnostic plots...\n")

# ── 1. Power curve ────────────────────────────────────────────────────────────
p_power <- ggplot(oc_df, aes(x = n_per_arm, y = Power,
                               colour = as.factor(dose),
                               fill   = as.factor(dose))) +
  geom_ribbon(aes(ymin = Power_CI_lo, ymax = Power_CI_hi), alpha = 0.15,
              colour = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.80, linetype = "dashed", colour = "tomato",
             linewidth = 0.9) +
  scale_colour_viridis_d(name = "Dose (mg)", option = "D") +
  scale_fill_viridis_d(name = "Dose (mg)",   option = "D") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_x_continuous(breaks = unique(oc_df$n_per_arm)) +
  labs(
    x       = "Sample size per arm (n)",
    y       = "Power",
    title   = "Power Curve — EC₅₀ Estimation Precision",
    caption = "Ribbon = Monte Carlo 95% CI (n = 1000 replicates). Dashed = 80% target."
  ) +
  theme(legend.position = "bottom")

save_figure(p_power, results_dir, "power_curve", dpi = 300, width = 9, height = 5)
cat("[03_diagnostics.R] Power curve saved.\n")

# ── 2. EC50 precision (boxplots) ──────────────────────────────────────────────
ec50_df <- raw_df %>%
  filter(dose > 0 & !is.na(EC50_est)) %>%
  mutate(n_label = paste0("n=", n_per_arm))

EC50_true <- 2.5

p_ec50 <- ggplot(ec50_df, aes(x = as.factor(n_per_arm), y = EC50_est,
                               fill = as.factor(dose))) +
  geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.3, position = position_dodge(0.8)) +
  geom_hline(yintercept = EC50_true, linetype = "dashed", colour = "tomato",
             linewidth = 0.9) +
  geom_hline(yintercept = EC50_true * 0.70, linetype = "dotted", colour = "grey50") +
  geom_hline(yintercept = EC50_true * 1.30, linetype = "dotted", colour = "grey50") +
  scale_fill_viridis_d(name = "Dose (mg)", option = "D") +
  scale_y_continuous(limits = c(0, pmin(quantile(ec50_df$EC50_est, 0.99, na.rm=TRUE), 15))) +
  labs(
    x       = "Sample size per arm",
    y       = "Estimated EC₅₀ (mg·h/L)",
    title   = "EC₅₀ Estimation Precision by Sample Size and Dose",
    caption = paste0("Dashed = true EC₅₀ (", EC50_true, "). ",
                     "Dotted = ±30% bounds. Points = outliers.")
  ) +
  theme(legend.position = "bottom")

save_figure(p_ec50, results_dir, "ec50_precision", dpi = 300, width = 10, height = 6)
cat("[03_diagnostics.R] EC50 precision plot saved.\n")

# ── 3. Operating characteristics table (gt) ───────────────────────────────────
oc_gt <- oc_df %>%
  mutate(
    Power_fmt   = sprintf("%.1f%% (%.1f–%.1f%%)",
                          Power*100, Power_CI_lo*100, Power_CI_hi*100),
    Bias_fmt    = sprintf("%+.1f%%", EC50_Bias * 100),
    RMSE_fmt    = sprintf("%.3f", EC50_RMSE),
    Cover_fmt   = sprintf("%.0f%%", EC50_Coverage * 100)
  ) %>%
  select(n_per_arm, dose, Power_fmt, Bias_fmt, RMSE_fmt, Cover_fmt) %>%
  rename(`n / arm` = n_per_arm, `Dose (mg)` = dose,
         `Power (95% CI)` = Power_fmt,
         `EC₅₀ Bias` = Bias_fmt,
         `EC₅₀ RMSE` = RMSE_fmt,
         `±30% Coverage` = Cover_fmt) |>
  gt(groupname_col = "Dose (mg)") |>
  tab_header(
    title    = "Operating Characteristics",
    subtitle = "Power, EC₅₀ bias, RMSE, and coverage by sample size and dose (n = 1000 replicates)"
  ) |>
  tab_style(
    style     = cell_fill(color = "#d4edda"),
    locations = cells_body(
      rows = as.numeric(gsub("%.*", "", `Power (95% CI)`)) >= 80
    )
  ) |>
  tab_footnote(
    "Green = ≥ 80% power. Coverage = % replicates with EC₅₀ within ±30% of true value."
  )

saveRDS(oc_gt, file.path(results_dir, "oc_table.rds"))
gt::gtsave(oc_gt, file.path(results_dir, "operating_characteristics_table.html"))
cat("[03_diagnostics.R] OC table saved.\n")

# ── 4. Bias vs n plot ─────────────────────────────────────────────────────────
p_bias <- ggplot(oc_df, aes(x = n_per_arm, y = EC50_Bias * 100,
                              colour = as.factor(dose))) +
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey50") +
  geom_hline(yintercept = c(-30, 30), linetype = "dashed", colour = "grey70") +
  geom_line(linewidth = 1.0) +
  geom_point(size = 3) +
  scale_colour_viridis_d(name = "Dose (mg)", option = "D") +
  scale_x_continuous(breaks = unique(oc_df$n_per_arm)) +
  labs(
    x       = "Sample size per arm",
    y       = "EC₅₀ estimation bias (%)",
    title   = "EC₅₀ Estimation Bias vs Sample Size",
    caption = "Dashed lines = ±30% bounds. Positive = overestimation."
  ) +
  theme(legend.position = "bottom")

save_figure(p_bias, results_dir, "ec50_bias", dpi = 300, width = 9, height = 5)
cat("[03_diagnostics.R] Bias plot saved.\n")

cat("[03_diagnostics.R] Project 05 diagnostics complete.\n")
