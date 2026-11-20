library(ggplot2)
library(data.table)

sim1   <- readRDS("firth/sim1.rds")
sim2   <- readRDS("firth/sim2.rds")
params <- readRDS("firth/sim_params.rds")
list2env(params, envir = environment())

# ── Plotting ──────────────────────────────────────────────────────────────────

make_plots <- function(dt, x_var, x_label, x_breaks, true_beta, log_x = FALSE) {
  # dt <- sim1
  # x_var <- "n"
  # x_label <- ""
  # x_breaks <- n_sizes
  # log_x = TRUE
  x_scale     <- if (log_x) scale_x_log10(breaks = x_breaks) else scale_x_continuous(breaks = x_breaks)
  
  # Plot 1: separation (data) + non-convergence (MLE and both Firth implementations)
  rates <- dt[, .(
    `Quasi-complete separation` = mean(quasi_complete,  na.rm = TRUE),
    `Near separation`           = mean(near_separated,  na.rm = TRUE),
    `MLE non-convergence`       = mean(mle_nonconv,     na.rm = TRUE),
    `Firth (brglm2) non-convergence`  = mean(firth_nonconv,   na.rm = TRUE),
    `Firth (logistf) non-convergence` = mean(logistf_nonconv, na.rm = TRUE)
  ), by = x_var]
  setnames(rates, x_var, "x")

  rates_long <- melt(rates, id.vars = "x",
                     variable.name = "rate_type", value.name = "rate")
  rates_long[, rate_type := factor(rate_type, levels = c(
    "Quasi-complete separation", "Near separation", "MLE non-convergence",
    "Firth (brglm2) non-convergence", "Firth (logistf) non-convergence"))]

  # note: quasi-complete separation and MLE non-convergence coincide (see note in qmd)
  p1 <- ggplot(rates_long, aes(x = x, y = rate, colour = rate_type,
                               group = rate_type)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    x_scale +
    scale_colour_manual(values = c(
      "Quasi-complete separation"       = "#2c7bb6",
      "Near separation"                 = "#abd9e9",
      "MLE non-convergence"             = "#d7191c",
      "Firth (brglm2) non-convergence"  = "#1a9641",
      "Firth (logistf) non-convergence" = "#fdae61"
    )) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
    labs(x = x_label, y = "% of simulations", colour = NULL) +
    theme_bw(base_size = 13) +
    theme(legend.position = "top", panel.grid = element_blank()) +
    guides(colour = guide_legend(nrow = 3))

  meth_lvls <- c("MLE", "Firth (brglm2)", "Firth (logistf)")
  meth_cols <- c("MLE" = "#E07B54", "Firth (brglm2)" = "#5B8DB8", "Firth (logistf)" = "#8FBF7F")

  # Plot 2: bias -- each method restricted to its own converged runs
  est_long <- rbindlist(list(
    dt[mle_nonconv == FALSE,     .(x = get(x_var), estimate = mle_est,     method = "MLE")],
    dt[firth_nonconv == FALSE,   .(x = get(x_var), estimate = firth_est,   method = "Firth (brglm2)")],
    dt[logistf_nonconv == FALSE, .(x = get(x_var), estimate = logistf_est, method = "Firth (logistf)")]
  ))
  est_long[, method := factor(method, levels = meth_lvls)]
  est_long[, bias   := estimate - true_beta]

  p2 <- ggplot(est_long, aes(x = x, y = bias, fill = method,
                              group = interaction(x, method))) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_boxplot() +
    x_scale +
    scale_fill_manual(values = meth_cols) +
    labs(
      x       = x_label,
      y       = expression("Bias = "~hat(beta) - beta),
      fill    = "Method"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "top", panel.grid = element_blank())

  # Plot 3: SE -- each method restricted to its own converged runs
  se_long <- rbindlist(list(
    dt[mle_nonconv == FALSE,     .(x = get(x_var), se = mle_se,     method = "MLE")],
    dt[firth_nonconv == FALSE,   .(x = get(x_var), se = firth_se,   method = "Firth (brglm2)")],
    dt[logistf_nonconv == FALSE, .(x = get(x_var), se = logistf_se, method = "Firth (logistf)")]
  ))
  se_long[, method := factor(method, levels = meth_lvls)]

  p3 <- ggplot(se_long, aes(x = x, y = se, fill = method,
                             group = interaction(x, method))) +
    geom_boxplot() +
    x_scale +
    scale_fill_manual(values = meth_cols) +
    labs(
      x       = x_label,
      y       = expression("SE("~hat(beta)~")"),
      fill    = "Method"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "top", panel.grid = element_blank())

  list(separation = p1, bias = p2, se = p3)
}

plots1 <- make_plots(sim1, "n",           "Sample size (n)",        n_sizes,           true_beta, log_x = TRUE)
plots2 <- make_plots(sim2, "prevalence",  "Overall prevalence (%)", prevalences * 100, true_beta, log_x = TRUE)
