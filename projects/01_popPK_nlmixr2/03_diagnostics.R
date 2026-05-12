# 03_diagnostics.R — Project 01: GOF, VPC, bootstrap, shrinkage, dosing plot
# Sources shared helpers from R/helpers.R. Reads fit object from results/.
# All figures saved as PNG (300 dpi) and SVG.

suppressPackageStartupMessages({
  library(nlmixr2)
  library(rxode2)
  library(dplyr)
  library(ggplot2)
  library(gt)
  library(here)
})

source(here("R", "helpers.R"))

results_dir <- here("projects", "01_popPK_nlmixr2", "results")
data_dir    <- here("projects", "01_popPK_nlmixr2", "data")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load fit and data ─────────────────────────────────────────────────────────
fit_path <- file.path(results_dir, "fit_saem.rds")
if (!file.exists(fit_path)) {
  stop("[03_diagnostics.R] fit_saem.rds not found — run 02_fit_model.R first")
}
fit     <- readRDS(fit_path)
pk_data <- read.csv(file.path(data_dir, "pk_data.csv"), stringsAsFactors = FALSE)

cat("[03_diagnostics.R] Fit loaded. Generating diagnostics...\n")

# ── 1. GOF plots ─────────────────────────────────────────────────────────────
gofs <- gof_plots(fit, title = "Population PK — 2-cmt Allometric Model")
save_gof_plots(gofs, results_dir, "gof", dpi = 300)
# Also save combined panel under the name expected by index.qmd
ggplot2::ggsave(file.path(results_dir, "gof_panel.png"),
                gofs$panel, dpi = 300, width = 11, height = 9)
ggplot2::ggsave(file.path(results_dir, "gof_panel.svg"),
                gofs$panel, width = 11, height = 9)
cat("[03_diagnostics.R] GOF plots saved.\n")

# ── 2. Prediction-corrected VPC by weight tertile ─────────────────────────────
# Stratification by BW tertile answers the clinical question: does the model
# adequately predict exposure in light, mid, and heavy patients separately?
pk_obs <- pk_data[pk_data$EVID == 0 & !is.na(pk_data$DV) & pk_data$DV > 0, ]
pk_obs$WTTERT <- pk_obs$BW_tertile

vpc_plot <- tryCatch(
  pcvpc(fit, pk_obs, strat_var = "WTTERT", n_sim = 500, seed = 42),
  error = function(e) {
    warning("pcvpc failed: ", conditionMessage(e), " — creating placeholder plot")
    ggplot(pk_obs, aes(x = TIME, y = DV, colour = BW_tertile)) +
      geom_point(alpha = 0.5) +
      facet_wrap(~BW_tertile) +
      scale_y_log10() +
      labs(title = "Observed concentrations by weight tertile (VPC pending)",
           x = "Time (h)", y = "Concentration (ng/mL)") +
      theme(legend.position = "none")
  }
)

save_figure(vpc_plot, results_dir, "pcvpc_by_weight", dpi = 300, width = 11, height = 5)
cat("[03_diagnostics.R] VPC saved.\n")

# ── 3. Bootstrap CI ───────────────────────────────────────────────────────────
n_boot  <- as.integer(Sys.getenv("BOOTSTRAP_N_BOOT", unset = "200"))
n_cores <- max(1L, parallel::detectCores() - 1L)

cat(sprintf("[03_diagnostics.R] Bootstrap: n=%d, cores=%d ...\n", n_boot, n_cores))

boot_df <- tryCatch(
  bootstrap_nlmixr(fit, pk_data, n_boot = n_boot, n_cores = n_cores, seed = 123),
  error = function(e) {
    warning("Bootstrap failed: ", conditionMessage(e))
    NULL
  }
)

if (!is.null(boot_df)) {
  boot_summary <- summarize_bootstrap(boot_df, ci = 0.95)
  saveRDS(boot_summary, file.path(results_dir, "bootstrap_summary.rds"))
  gt::gtsave(boot_summary, file.path(results_dir, "bootstrap_summary.html"))
  write.csv(as.data.frame(boot_df), file.path(results_dir, "bootstrap_raw.csv"),
            row.names = FALSE)
  cat("[03_diagnostics.R] Bootstrap CI table saved.\n")
}

# ── 4. Shrinkage table ────────────────────────────────────────────────────────
shrink_tbl <- tryCatch(
  shrinkage_table(fit, threshold_eta = 0.30, threshold_eps = 0.20),
  error = function(e) {
    warning("shrinkage_table failed: ", conditionMessage(e))
    NULL
  }
)

if (!is.null(shrink_tbl)) {
  saveRDS(shrink_tbl, file.path(results_dir, "shrinkage_table.rds"))
  gt::gtsave(shrink_tbl, file.path(results_dir, "shrinkage_table.html"))
  cat("[03_diagnostics.R] Shrinkage table saved.\n")
}

# ── 5. Condition number ───────────────────────────────────────────────────────
# Condition number of the variance-covariance matrix: sqrt(lambda_max / lambda_min).
# Values > 1000 suggest near-singular covariance, indicating parameter correlations
# that may make individual parameters poorly identifiable.
cond_num <- tryCatch({
  vcv        <- fit$cov
  eigenvals  <- eigen(vcv)$values
  eigenvals  <- eigenvals[eigenvals > 0]
  sqrt(max(eigenvals) / min(eigenvals))
}, error = function(e) NA_real_)

cat(sprintf("[03_diagnostics.R] Condition number: %.1f%s\n",
            cond_num, if (!is.na(cond_num) && cond_num > 1000) " [HIGH — check correlations]" else ""))

writeLines(sprintf("Condition number: %.2f", cond_num),
           file.path(results_dir, "condition_number.txt"))

# ── 6. Dosing comparison plot ─────────────────────────────────────────────────
dose_sim_path <- file.path(results_dir, "dosing_simulation.csv")
if (file.exists(dose_sim_path)) {
  dose_sim <- read.csv(dose_sim_path)
  TARGET   <- 1.0  # ng/mL

  # Reshape to long format for faceted plot
  dose_long <- tidyr::pivot_longer(
    dose_sim,
    cols      = c(Trough_flat, Trough_banded),
    names_to  = "Strategy",
    values_to = "Trough_ngmL"
  ) %>%
    mutate(Strategy = recode(Strategy,
                              Trough_flat   = "Flat 100 mg",
                              Trough_banded = "Weight-banded\n(75/100/125 mg)"))

  p_dose <- ggplot(dose_long,
                   aes(x = BW, y = Trough_ngmL, colour = BW_group)) +
    geom_line(linewidth = 1.0) +
    geom_hline(yintercept = TARGET, linetype = "dashed", colour = "tomato",
               linewidth = 0.9) +
    facet_wrap(~Strategy) +
    scale_colour_viridis_d(name = "Weight group") +
    labs(x = "Body weight (kg)",
         y = "Predicted steady-state trough (ng/mL)",
         title = "Flat vs Weight-banded Dosing — Simulated Trough Exposures",
         caption = "Dashed line = target trough (1 ng/mL, illustrative)") +
    theme(legend.position = "bottom")

  save_figure(p_dose, results_dir, "dosing_comparison", dpi = 300, width = 10, height = 5)
  cat("[03_diagnostics.R] Dosing comparison plot saved.\n")
}

cat("[03_diagnostics.R] All diagnostics complete.\n")
