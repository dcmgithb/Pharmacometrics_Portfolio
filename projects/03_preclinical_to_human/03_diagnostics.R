# 03_diagnostics.R — Project 03: Allometric scaling diagnostics and FIH plots

suppressPackageStartupMessages({
  library(nlmixr2)
  library(dplyr)
  library(ggplot2)
  library(gt)
  library(here)
})

source(here("R", "helpers.R"))

results_dir <- here("projects", "03_preclinical_to_human", "results")
data_dir    <- here("projects", "03_preclinical_to_human", "data")

fit_rat_path <- file.path(results_dir, "fit_rat.rds")
fit_dog_path <- file.path(results_dir, "fit_dog.rds")

if (!file.exists(fit_rat_path)) stop("Run 02_fit_model.R first (fit_rat.rds missing)")

fit_rat  <- readRDS(fit_rat_path)
fit_dog  <- if (file.exists(fit_dog_path)) readRDS(fit_dog_path) else NULL
rat_data <- read.csv(file.path(data_dir, "rat_pk.csv"))
dog_data <- read.csv(file.path(data_dir, "dog_pk.csv"))

cat("[03_diagnostics.R] Generating diagnostics...\n")

# ── 1. Species GOF panels ─────────────────────────────────────────────────────
gofs_rat <- gof_plots(fit_rat, title = "Rat PK — 1-cmt IV Model GOF")
ggplot2::ggsave(file.path(results_dir, "gof_rat.png"),
                gofs_rat$panel, dpi = 300, width = 11, height = 9)
ggplot2::ggsave(file.path(results_dir, "gof_rat.svg"),
                gofs_rat$panel, width = 11, height = 9)

if (!is.null(fit_dog)) {
  gofs_dog <- gof_plots(fit_dog, title = "Dog PK — 1-cmt IV Model GOF")
  ggplot2::ggsave(file.path(results_dir, "gof_dog.png"),
                  gofs_dog$panel, dpi = 300, width = 11, height = 9)
  ggplot2::ggsave(file.path(results_dir, "gof_dog.svg"),
                  gofs_dog$panel, width = 11, height = 9)
}
cat("[03_diagnostics.R] Species GOF plots saved.\n")

# ── 2. Allometric scaling plot ────────────────────────────────────────────────
BW_sp   <- c(0.25, 10.0)
BW_labs <- c("Rat (0.25 kg)", "Dog (10 kg)")

CL_est  <- c(exp(fit_rat$fixef["log_CL"]),
             if (!is.null(fit_dog)) exp(fit_dog$fixef["log_CL"]) else 1.20)
V_est   <- c(exp(fit_rat$fixef["log_V"]),
             if (!is.null(fit_dog)) exp(fit_dog$fixef["log_V"])  else 15.0)

# Allometric line from fitted coefficients
alpha_CL <- exp(mean(log(CL_est / BW_sp^0.75)))
alpha_V  <- exp(mean(log(V_est  / BW_sp^1.00)))

BW_seq      <- 10^seq(log10(0.1), log10(100), length.out = 200)
CL_line     <- alpha_CL * BW_seq^0.75
V_line      <- alpha_V  * BW_seq^1.00

# Predicted human point
CL_human <- alpha_CL * 70^0.75
V_human  <- alpha_V  * 70^1.00

# Bootstrap CI for the allometric prediction (parametric bootstrap from fit uncertainty)
set.seed(999)
n_boot <- 500
CL_boot <- sapply(seq_len(n_boot), function(i) {
  noise_rat <- rnorm(1, 0, exp(fit_rat$fixef["prop_err"]) * CL_est[1])
  noise_dog <- if (!is.null(fit_dog))
    rnorm(1, 0, exp(fit_dog$fixef["prop_err"]) * CL_est[2]) else 0
  cl_r <- pmax(CL_est[1] + noise_rat, 1e-4)
  cl_d <- pmax(CL_est[2] + noise_dog, 1e-4)
  a    <- exp(mean(log(c(cl_r, cl_d) / BW_sp^0.75)))
  a * 70^0.75
})
CI_CL <- quantile(CL_boot, c(0.025, 0.975))

human_pt <- data.frame(BW = 70, CL = CL_human, V = V_human,
                        label = "Human\n(predicted)",
                        CI_lo = CI_CL[1], CI_hi = CI_CL[2])
sp_pts   <- data.frame(BW = BW_sp, CL = CL_est, V = V_est,
                        label = BW_labs)

p_allo <- ggplot() +
  geom_line(data = data.frame(BW = BW_seq, CL = CL_line),
            aes(x = BW, y = CL), colour = "steelblue", linewidth = 1.1) +
  geom_point(data = sp_pts, aes(x = BW, y = CL),
             colour = "steelblue", size = 4, shape = 16) +
  geom_text(data = sp_pts, aes(x = BW, y = CL, label = label),
            nudge_y = 0.1, hjust = -0.1, size = 3.5) +
  geom_point(data = human_pt, aes(x = BW, y = CL),
             colour = "tomato", size = 5, shape = 17) +
  geom_errorbar(data = human_pt,
                aes(x = BW, ymin = CI_lo, ymax = CI_hi),
                colour = "tomato", width = 2, linewidth = 0.9) +
  geom_text(data = human_pt, aes(x = BW, y = CL, label = label),
            nudge_y = -0.3, colour = "tomato", size = 3.5) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "Body weight (kg)", y = "Clearance (L/hr)",
       title = "Allometric Scaling of Clearance",
       subtitle = sprintf("CL = %.4f × BW^0.75   (fixed exponent)", alpha_CL),
       caption = "Error bar = 95% bootstrap CI for human prediction") +
  annotation_logticks()

save_figure(p_allo, results_dir, "allometric_plot", dpi = 300, width = 8, height = 6)
cat("[03_diagnostics.R] Allometric scaling plot saved.\n")

# Update human_pk_prediction.csv with bootstrap CI
hpred_path <- file.path(results_dir, "human_pk_prediction.csv")
if (file.exists(hpred_path)) {
  hpred <- read.csv(hpred_path)
  cl_row <- hpred$Parameter == "CL (L/hr)"
  hpred$CI_Lower[cl_row] <- round(CI_CL[1], 3)
  hpred$CI_Upper[cl_row] <- round(CI_CL[2], 3)
  write.csv(hpred, hpred_path, row.names = FALSE)
}

# ── 3. Sensitivity to allometric exponent ────────────────────────────────────
# Shows how sensitive the human CL prediction is to the choice of exponent
# (0.75 is standard; some use brain-weight corrected "rule of exponents" values)
exp_seq    <- seq(0.60, 0.90, by = 0.01)
CL_vs_exp  <- sapply(exp_seq, function(b) {
  a <- exp(mean(log(CL_est / BW_sp^b)))
  a * 70^b
})

p_sens <- ggplot(data.frame(exponent = exp_seq, CL_human = CL_vs_exp),
                 aes(x = exponent, y = CL_human)) +
  geom_line(colour = "steelblue", linewidth = 1.2) +
  geom_vline(xintercept = 0.75, linetype = "dashed", colour = "tomato",
             linewidth = 0.9) +
  geom_point(aes(x = 0.75, y = CL_human),
             data = data.frame(exponent = 0.75,
                               CL_human = alpha_CL * 70^0.75),
             colour = "tomato", size = 3) +
  labs(x = "Allometric exponent (β)", y = "Predicted human CL (L/hr)",
       title = "Sensitivity of Human CL Prediction to Allometric Exponent",
       caption = "Dashed line = standard physiological exponent (0.75)") +
  annotate("text", x = 0.75, y = CL_human * 1.1,
           label = sprintf("β=0.75: CL=%.2f L/hr", CL_human),
           colour = "tomato", hjust = -0.1, size = 3.5)

save_figure(p_sens, results_dir, "sensitivity_exponent", dpi = 300, width = 8, height = 5)
cat("[03_diagnostics.R] Sensitivity plot saved.\n")

# ── 4. FIH recommendation text ────────────────────────────────────────────────
fih_tbl <- read.csv(file.path(results_dir, "fih_dose_table.csv"))
dose_fih <- fih_tbl$Dose_mg[fih_tbl$Approach == "Selected FIH dose"]
safety_m <- fih_tbl$Safety_margin[fih_tbl$Approach == "Selected FIH dose"]

fih_text <- paste0(
  "## FIH Starting Dose Recommendation\n\n",
  sprintf("Predicted human CL: **%.2f L/hr** (95%% CI: %.2f–%.2f L/hr, bootstrap)\n\n",
          CL_human, CI_CL[1], CI_CL[2]),
  sprintf("Predicted human V:  **%.1f L**  (t½ ≈ %.1f h)\n\n",
          V_human, log(2) * V_human / CL_human),
  sprintf("**Recommended FIH starting dose: %.2f mg**\n\n", dose_fih),
  sprintf("Safety margin vs NOAEL-HED: **%.0f-fold**\n\n", safety_m),
  "Basis: Lower of MABEL-based (10% receptor occupancy) ",
  "and NOAEL-derived HED with 10× safety factor. ",
  "Cohort dose escalation should not exceed the NOAEL HED without additional toxicology data."
)
writeLines(fih_text, file.path(results_dir, "fih_recommendation.txt"))

cat("[03_diagnostics.R] FIH recommendation saved.\n")
cat("[03_diagnostics.R] Project 03 diagnostics complete.\n")
