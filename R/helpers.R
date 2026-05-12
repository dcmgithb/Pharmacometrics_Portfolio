# R/helpers.R — shared pharmacometric diagnostic functions
# Sourced by every project's 03_diagnostics.R script.
#
# Functions exported:
#   gof_plots()           — standard GOF panel (4 scatter plots)
#   save_gof_plots()      — save GOF list to PNG + SVG
#   pcvpc()               — prediction-corrected VPC with stratification
#   bootstrap_nlmixr()    — case-resampling bootstrap CI
#   summarize_bootstrap() — bootstrap results → gt table
#   shrinkage_table()     — η- and ε-shrinkage → gt table
#   save_figure()         — save any ggplot as PNG + SVG
#   set_portfolio_theme() — apply consistent ggplot2 theme globally

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(scales)
  library(gt)
  library(future)
  library(future.apply)
  library(vpc)
})

# ── Theme ─────────────────────────────────────────────────────────────────────

set_portfolio_theme <- function() {
  theme_set(
    theme_bw(base_size = 12) +
      theme(
        strip.background = element_rect(fill = "grey92"),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom"
      )
  )
  # Colorblind-safe palette for discrete scales
  options(ggplot2.discrete.colour   = scales::hue_pal()(8),
          ggplot2.discrete.fill     = scales::hue_pal()(8))
  invisible(NULL)
}

set_portfolio_theme()

# ── Figure saving ─────────────────────────────────────────────────────────────

#' Save a ggplot as both PNG and SVG
#'
#' @param plot    ggplot object
#' @param dir     output directory (created if absent)
#' @param name    base filename without extension
#' @param dpi     PNG resolution; SVG is vector so dpi ignored
#' @param width   inches
#' @param height  inches
save_figure <- function(plot, dir, name, dpi = 300, width = 8, height = 6) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  png_path <- file.path(dir, paste0(name, ".png"))
  svg_path <- file.path(dir, paste0(name, ".svg"))
  ggplot2::ggsave(png_path, plot, dpi = dpi, width = width, height = height)
  ggplot2::ggsave(svg_path, plot, width = width, height = height)
  invisible(list(png = png_path, svg = svg_path))
}

# ── GOF plots ─────────────────────────────────────────────────────────────────

#' Standard goodness-of-fit diagnostic panel
#'
#' @param fit       nlmixr2 fit object OR data.frame with DV, PRED, IPRED,
#'                  CWRES, TIME columns.
#' @param title     Panel title string.
#' @param dv_var    Column name for observed DV.
#' @param pred_var  Column name for population prediction.
#' @param ipred_var Column name for individual prediction.
#' @param time_var  Column name for time.
#' @param cwres_var Column name for CWRES.
#' @param log_scale TRUE, FALSE, or "auto". "auto" uses log scale when
#'                  DV spans > 2 orders of magnitude — prevents visual
#'                  compression at low concentrations.
#'
#' @return Named list: dv_pred, dv_ipred, cwres_time, cwres_pred, panel.
gof_plots <- function(fit,
                      title     = "GOF Diagnostics",
                      dv_var    = "DV",
                      pred_var  = "PRED",
                      ipred_var = "IPRED",
                      time_var  = "TIME",
                      cwres_var = "CWRES",
                      log_scale = "auto") {

  # Accept either a nlmixr2 fit or a plain data.frame
  if (inherits(fit, "nlmixr2FitCore")) {
    df <- as.data.frame(fit)
  } else {
    df <- fit
  }

  required <- c(dv_var, pred_var, ipred_var, time_var, cwres_var)
  missing  <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("gof_plots: missing columns: ", paste(missing, collapse = ", "))
  }

  obs <- df[df[[dv_var]] > 0 & !is.na(df[[dv_var]]) & !is.na(df[[cwres_var]]), ]

  # Decide whether to use log scale
  use_log <- if (isTRUE(log_scale == "auto")) {
    range_ratio <- max(obs[[dv_var]], na.rm = TRUE) /
                   max(min(obs[[dv_var]][obs[[dv_var]] > 0], na.rm = TRUE), 1e-9)
    range_ratio > 100
  } else {
    isTRUE(log_scale)
  }

  scale_layer <- if (use_log) {
    list(scale_x_log10(), scale_y_log10())
  } else {
    list()
  }

  # DV vs PRED
  p_dv_pred <- ggplot(obs, aes(x = .data[[pred_var]], y = .data[[dv_var]],
                               colour = abs(.data[[cwres_var]]))) +
    geom_abline(slope = 1, intercept = 0, colour = "tomato", linewidth = 0.8) +
    geom_point(alpha = 0.45, size = 1.6) +
    geom_smooth(method = "loess", se = FALSE, colour = "steelblue",
                linewidth = 0.8, formula = y ~ x) +
    scale_colour_viridis_c(name = "|CWRES|", option = "D", direction = -1) +
    scale_layer +
    labs(x = "Population prediction (PRED)", y = "Observed (DV)",
         title = "DV vs PRED") +
    theme(legend.position = "right")

  # DV vs IPRED
  p_dv_ipred <- ggplot(obs, aes(x = .data[[ipred_var]], y = .data[[dv_var]],
                                colour = abs(.data[[cwres_var]]))) +
    geom_abline(slope = 1, intercept = 0, colour = "tomato", linewidth = 0.8) +
    geom_point(alpha = 0.45, size = 1.6) +
    geom_smooth(method = "loess", se = FALSE, colour = "steelblue",
                linewidth = 0.8, formula = y ~ x) +
    scale_colour_viridis_c(name = "|CWRES|", option = "D", direction = -1) +
    scale_layer +
    labs(x = "Individual prediction (IPRED)", y = "Observed (DV)",
         title = "DV vs IPRED") +
    theme(legend.position = "right")

  # Reference lines for CWRES plots: ±1.96 because under a correctly specified
  # model CWRES ~ N(0,1), so 95% of values should fall within ±1.96.
  cwres_refs <- list(
    geom_hline(yintercept = 0,    linetype = "solid",  colour = "grey50"),
    geom_hline(yintercept =  1.96, linetype = "dashed", colour = "tomato", alpha = 0.7),
    geom_hline(yintercept = -1.96, linetype = "dashed", colour = "tomato", alpha = 0.7)
  )

  # CWRES vs TIME
  p_cwres_time <- ggplot(obs, aes(x = .data[[time_var]], y = .data[[cwres_var]])) +
    cwres_refs +
    geom_point(alpha = 0.4, size = 1.4, colour = "steelblue") +
    geom_smooth(method = "loess", se = FALSE, colour = "darkorange",
                linewidth = 0.8, formula = y ~ x) +
    coord_cartesian(ylim = c(min(-4, min(obs[[cwres_var]], na.rm=TRUE)),
                             max( 4, max(obs[[cwres_var]], na.rm=TRUE)))) +
    labs(x = "Time (h)", y = "CWRES", title = "CWRES vs TIME")

  # CWRES vs PRED
  p_cwres_pred <- ggplot(obs, aes(x = .data[[pred_var]], y = .data[[cwres_var]])) +
    cwres_refs +
    geom_point(alpha = 0.4, size = 1.4, colour = "steelblue") +
    geom_smooth(method = "loess", se = FALSE, colour = "darkorange",
                linewidth = 0.8, formula = y ~ x) +
    coord_cartesian(ylim = c(min(-4, min(obs[[cwres_var]], na.rm=TRUE)),
                             max( 4, max(obs[[cwres_var]], na.rm=TRUE)))) +
    labs(x = "Population prediction (PRED)", y = "CWRES", title = "CWRES vs PRED")

  panel <- cowplot::plot_grid(p_dv_pred, p_dv_ipred, p_cwres_time, p_cwres_pred,
                              ncol = 2, labels = "AUTO")
  panel <- cowplot::ggdraw(panel) +
    cowplot::draw_label(title, x = 0.5, y = 1.01, hjust = 0.5, vjust = 0,
                        fontface = "bold", size = 13)

  list(dv_pred    = p_dv_pred,
       dv_ipred   = p_dv_ipred,
       cwres_time = p_cwres_time,
       cwres_pred = p_cwres_pred,
       panel      = panel)
}

#' Save GOF plot list to PNG and SVG
#'
#' @param gof_list List returned by gof_plots().
#' @param dir      Output directory.
#' @param prefix   Filename prefix.
#' @param dpi      PNG resolution.
save_gof_plots <- function(gof_list, dir, prefix, dpi = 300) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  nms <- c("dv_pred", "dv_ipred", "cwres_time", "cwres_pred")
  for (nm in nms) {
    save_figure(gof_list[[nm]], dir, paste0(prefix, "_", nm), dpi = dpi,
                width = 6, height = 5)
  }
  # Panel: wider + taller
  save_figure(gof_list$panel, dir, paste0(prefix, "_panel"), dpi = dpi,
              width = 11, height = 9)
  invisible(NULL)
}

# ── Prediction-corrected VPC ──────────────────────────────────────────────────

#' Generate a prediction-corrected VPC
#'
#' @param fit           nlmixr2 fit object.
#' @param data          Original dataset used for fitting.
#' @param strat_var     Column name for stratification facets, or NULL.
#' @param n_sim         Number of simulated datasets.
#' @param ci            Confidence interval width for ribbon (e.g. 0.95).
#' @param obs_pi        Observed percentiles to overlay (length-3 vector).
#' @param seed          RNG seed for simulation.
#' @param n_bins        Number of time bins for binning.
#'
#' @details
#' Prediction-corrected VPC divides each concentration by the median
#' population prediction in its time bin. This removes the structural model's
#' covariate effects and allows pooling across dose/covariate groups without
#' bias — something a standard VPC cannot do cleanly.
#'
#' @return ggplot object.
pcvpc <- function(fit,
                  data,
                  strat_var = NULL,
                  n_sim     = 1000,
                  ci        = 0.95,
                  obs_pi    = c(0.05, 0.50, 0.95),
                  seed      = 42,
                  n_bins    = 10) {

  set.seed(seed)

  # Extract fit data if nlmixr2 object passed
  fit_df <- if (inherits(fit, "nlmixr2FitCore")) as.data.frame(fit) else fit

  obs <- fit_df[fit_df$EVID == 0 & fit_df$DV > 0 & !is.na(fit_df$DV), ]

  # Simulate from model using rxode2 if available; else use vpc package directly
  if (requireNamespace("vpc", quietly = TRUE)) {
    strat_col <- if (!is.null(strat_var) && strat_var %in% names(obs)) strat_var else NULL

    # Build simulated dataset by re-using the vpc package's simulation wrapper
    # vpc_vpc returns a VPC-ready data structure
    tryCatch({
      sim_data <- vpc::sim_data(
        m    = fit,
        nsim = n_sim,
        dv   = "DV",
        obs_col = "DV",
        id_col  = "ID",
        idv_col = "TIME"
      )
    }, error = function(e) {
      # Fall back to parametric simulation from model parameters
      sim_data <- NULL
    })

    if (is.null(sim_data)) {
      # Minimal fallback: use observed data as both obs and sim (displays structure)
      warning("pcvpc: could not simulate from model — displaying observed data only")
      p <- ggplot(obs, aes(x = TIME, y = DV)) +
        geom_point(alpha = 0.4, colour = "steelblue") +
        labs(x = "Time (h)", y = "Concentration",
             title = "VPC (simulation unavailable — run 02_fit_model.R first)",
             caption = "Dashed lines = observed 5th/50th/95th percentiles") +
        stat_summary(fun = median, geom = "line", colour = "darkred", linewidth = 0.8)
      if (!is.null(strat_col)) p <- p + facet_wrap(as.formula(paste0("~", strat_col)))
      return(p)
    }

    vpc_plot <- vpc::vpc(sim = sim_data, obs = obs,
                         obs_cols  = list(dv = "DV", idv = "TIME", id = "ID"),
                         sim_cols  = list(dv = "DV", idv = "TIME", id = "ID"),
                         pi        = obs_pi[c(1,3)],
                         ci        = c((1-ci)/2, 1-(1-ci)/2),
                         pred_corr = TRUE,
                         bins      = n_bins,
                         stratify  = strat_col,
                         show      = list(obs_dv = TRUE, obs_ci = TRUE,
                                          pi = TRUE, pi_ci = TRUE))
    return(vpc_plot)
  }

  # Plain percentile overlay if vpc package unavailable
  alpha_ci <- (1 - ci) / 2
  time_breaks <- quantile(obs$TIME, probs = seq(0, 1, length.out = n_bins + 1))
  obs$bin <- cut(obs$TIME, breaks = unique(time_breaks), include.lowest = TRUE)

  pctiles <- obs |>
    group_by(bin) |>
    summarise(
      t_mid = median(TIME),
      p05   = quantile(DV, obs_pi[1], na.rm = TRUE),
      p50   = quantile(DV, obs_pi[2], na.rm = TRUE),
      p95   = quantile(DV, obs_pi[3], na.rm = TRUE),
      .groups = "drop"
    )

  p <- ggplot(obs, aes(x = TIME, y = DV)) +
    geom_point(alpha = 0.3, size = 1.2, colour = "grey60") +
    geom_line(data = pctiles, aes(x = t_mid, y = p05), linetype = "dashed",
              colour = "steelblue") +
    geom_line(data = pctiles, aes(x = t_mid, y = p50), colour = "darkred", linewidth = 1) +
    geom_line(data = pctiles, aes(x = t_mid, y = p95), linetype = "dashed",
              colour = "steelblue") +
    labs(x = "Time (h)", y = "Concentration",
         title = "Observed percentile overlay (vpc package not available)",
         caption = "Dashed = 5th/95th percentile; solid = median")

  if (!is.null(strat_var) && strat_var %in% names(obs)) {
    p <- p + facet_wrap(as.formula(paste0("~", strat_var)))
  }
  p
}

# ── Bootstrap ─────────────────────────────────────────────────────────────────

#' Case-resampling bootstrap CI for nlmixr2 parameter estimates
#'
#' @param fit      nlmixr2 fit object (converged).
#' @param data     Dataset used for fitting; must have an ID column.
#' @param n_boot   Number of bootstrap replicates.
#' @param n_cores  Number of parallel workers.
#' @param seed     Base RNG seed (replicate i uses seed + i).
#' @param id_var   Name of the subject identifier column.
#'
#' @details
#' Case resampling (resample subjects with replacement) makes no distributional
#' assumption about residuals, making it robust to model misspecification.
#' This is the standard approach per EMA (2007) and FDA (2019) PopPK guidelines.
#'
#' @return Data.frame with n_boot rows, one column per parameter plus
#'         boot_id and convergence_flag.
bootstrap_nlmixr <- function(fit,
                             data,
                             n_boot   = 200,
                             n_cores  = max(1L, parallel::detectCores() - 1L),
                             seed     = 123,
                             id_var   = "ID") {

  n_boot_env <- as.integer(Sys.getenv("BOOTSTRAP_N_BOOT", unset = as.character(n_boot)))
  if (n_boot_env != n_boot) {
    message(sprintf("bootstrap_nlmixr: BOOTSTRAP_N_BOOT=%d overrides n_boot=%d",
                    n_boot_env, n_boot))
    n_boot <- n_boot_env
  }

  ids     <- unique(data[[id_var]])
  n_subj  <- length(ids)

  # Extract model specification from the fit for re-fitting
  # Store original control settings
  orig_model   <- fit$model
  orig_control <- fit$control

  future::plan(future::multisession, workers = n_cores)
  on.exit(future::plan(future::sequential), add = TRUE)

  results <- future.apply::future_lapply(
    seq_len(n_boot),
    function(i) {
      set.seed(seed + i)
      sampled_ids <- sample(ids, n_subj, replace = TRUE)

      # Build resampled dataset: duplicate rows for each sampled ID,
      # reassigning a new unique ID to avoid duplicate-ID issues in nlmixr2
      boot_data <- do.call(rbind, lapply(seq_along(sampled_ids), function(j) {
        rows        <- data[data[[id_var]] == sampled_ids[j], ]
        rows[[id_var]] <- j  # new sequential ID
        rows
      }))

      fit_boot <- tryCatch({
        nlmixr2::nlmixr2(orig_model, boot_data,
                         est     = "saem",
                         control = orig_control)
      }, error = function(e) NULL)

      if (is.null(fit_boot)) {
        return(data.frame(boot_id = i, convergence_flag = 1L))
      }

      params <- as.data.frame(t(fit_boot$theta))
      params$boot_id          <- i
      params$convergence_flag <- 0L
      params
    },
    future.seed = TRUE
  )

  out <- do.call(dplyr::bind_rows, results)
  out <- out[order(out$boot_id), ]
  out
}

#' Summarise bootstrap results as a gt table
#'
#' @param boot_df  Data.frame from bootstrap_nlmixr().
#' @param orig_fit nlmixr2 fit object for point estimates.
#' @param ci       CI coverage (e.g. 0.95 for 95% CI).
#' @param exclude_nonconverged Logical; exclude rows with convergence_flag != 0.
#'
#' @return gt table object.
summarize_bootstrap <- function(boot_df,
                                orig_fit              = NULL,
                                ci                    = 0.95,
                                exclude_nonconverged  = TRUE) {

  if (exclude_nonconverged && "convergence_flag" %in% names(boot_df)) {
    n_fail <- sum(boot_df$convergence_flag != 0, na.rm = TRUE)
    if (n_fail > 0) message(sprintf("Excluding %d non-converged replicates", n_fail))
    boot_df <- boot_df[boot_df$convergence_flag == 0, ]
  }

  param_cols <- setdiff(names(boot_df), c("boot_id", "convergence_flag"))
  alpha      <- (1 - ci) / 2

  summ <- lapply(param_cols, function(p) {
    x <- boot_df[[p]]
    data.frame(
      Parameter     = p,
      Boot_Mean     = mean(x,                   na.rm = TRUE),
      Boot_Median   = median(x,                 na.rm = TRUE),
      CI_Lower      = quantile(x, alpha,         na.rm = TRUE),
      CI_Upper      = quantile(x, 1 - alpha,     na.rm = TRUE),
      Boot_SE       = sd(x,                      na.rm = TRUE),
      N_replicates  = sum(!is.na(x))
    )
  }) |> do.call(rbind)

  gt_tbl <- gt::gt(summ) |>
    gt::tab_header(
      title    = sprintf("Bootstrap Parameter Summary (%.0f%% CI)", ci * 100),
      subtitle = sprintf("n = %d replicates, case resampling", nrow(boot_df))
    ) |>
    gt::fmt_number(columns = c(Boot_Mean, Boot_Median, CI_Lower, CI_Upper, Boot_SE),
                   decimals = 4) |>
    gt::fmt_integer(columns = N_replicates)

  gt_tbl
}

# ── Shrinkage table ───────────────────────────────────────────────────────────

#' η- and ε-shrinkage table from a nlmixr2 fit
#'
#' @param fit           nlmixr2 fit object.
#' @param threshold_eta Shrinkage threshold for colour coding. Default 0.30.
#'                      Different regulatory thresholds (EMA: ~0.25, FDA: ~0.30)
#'                      can be set here to match the intended audience.
#' @param threshold_eps ε-shrinkage alert level. Default 0.20.
#'
#' @details
#' η-shrinkage = 1 - SD(η_i) / ω. High shrinkage (> threshold) means individual
#' empirical Bayes estimates (EBEs) collapse toward the population mean, making
#' them unreliable for covariate model building or individual dosing.
#'
#' ε-shrinkage = 1 - SD(IWRES). High ε-shrinkage implies sparse data; CWRES
#' distributions and GOF plots become uninformative.
#'
#' @return gt table object with colour-coded shrinkage column.
shrinkage_table <- function(fit,
                            threshold_eta = 0.30,
                            threshold_eps = 0.20) {

  if (!inherits(fit, "nlmixr2FitCore")) {
    stop("shrinkage_table: fit must be a nlmixr2FitCore object")
  }

  fit_df <- as.data.frame(fit)

  # η-shrinkage for each random effect
  eta_cols <- grep("^eta\\.", names(fit_df), value = TRUE)
  omega     <- sqrt(diag(fit$omega))  # population SD (sqrt of diagonal variances)

  eta_rows <- if (length(eta_cols) > 0 && length(omega) > 0) {
    n_etas <- min(length(eta_cols), length(omega))
    lapply(seq_len(n_etas), function(i) {
      sd_eta  <- sd(fit_df[[eta_cols[i]]], na.rm = TRUE)
      shrink  <- 1 - sd_eta / omega[i]
      data.frame(
        Parameter  = eta_cols[i],
        Type       = "η (IIV)",
        Omega_SD   = round(omega[i], 4),
        Shrinkage  = round(shrink, 4),
        stringsAsFactors = FALSE
      )
    }) |> do.call(rbind)
  } else {
    data.frame(Parameter = character(0), Type = character(0),
               Omega_SD = numeric(0), Shrinkage = numeric(0))
  }

  # ε-shrinkage from individual weighted residuals
  eps_shrink <- if ("IWRES" %in% names(fit_df)) {
    1 - sd(fit_df$IWRES[fit_df$EVID == 0], na.rm = TRUE)
  } else if ("CWRES" %in% names(fit_df)) {
    1 - sd(fit_df$CWRES[fit_df$EVID == 0], na.rm = TRUE)
  } else {
    NA_real_
  }

  eps_row <- data.frame(
    Parameter  = "epsilon (RUV)",
    Type       = "ε (residual)",
    Omega_SD   = NA_real_,
    Shrinkage  = round(eps_shrink, 4),
    stringsAsFactors = FALSE
  )

  tbl <- rbind(eta_rows, eps_row)

  # Interpretation labels
  tbl$Interpretation <- dplyr::case_when(
    is.na(tbl$Shrinkage)              ~ "N/A",
    tbl$Shrinkage < threshold_eta     ~ "Low — EBEs reliable",
    tbl$Shrinkage < 2 * threshold_eta ~ "Moderate — use with caution",
    TRUE                               ~ "High — EBEs unreliable"
  )

  tbl$Shrinkage_pct <- sprintf("%.1f%%", tbl$Shrinkage * 100)

  gt_tbl <- gt::gt(tbl[, c("Parameter", "Type", "Omega_SD", "Shrinkage_pct", "Interpretation")]) |>
    gt::tab_header(
      title    = "Shrinkage Summary",
      subtitle = sprintf("Alert thresholds: η > %.0f%%, ε > %.0f%%",
                         threshold_eta * 100, threshold_eps * 100)
    ) |>
    gt::cols_label(
      Omega_SD      = "ω or σ (SD)",
      Shrinkage_pct = "Shrinkage (%)",
      Interpretation = "Interpretation"
    ) |>
    gt::tab_style(
      style     = gt::cell_fill(color = "#d4edda"),
      locations = gt::cells_body(
        rows = grepl("Low", tbl$Interpretation)
      )
    ) |>
    gt::tab_style(
      style     = gt::cell_fill(color = "#fff3cd"),
      locations = gt::cells_body(
        rows = grepl("Moderate", tbl$Interpretation)
      )
    ) |>
    gt::tab_style(
      style     = gt::cell_fill(color = "#f8d7da"),
      locations = gt::cells_body(
        rows = grepl("High", tbl$Interpretation)
      )
    ) |>
    gt::fmt_number(columns = "Omega_SD", decimals = 4, rows = !is.na(tbl$Omega_SD)) |>
    gt::tab_footnote(
      footnote = paste0(
        "η-shrinkage = 1 − SD(η_i)/ω. ε-shrinkage = 1 − SD(IWRES). ",
        "High shrinkage implies EBEs and GOF plots are less informative. ",
        "Threshold: EMA guidance ~25%, FDA informal ~30%."
      )
    )

  gt_tbl
}
