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
  
  # Plot 1: separation (data) + non-convergence (MLE and Firth)
  rates <- dt[, .(
    `Complete separation`  = mean(separated,      na.rm = TRUE),
    `Near/quasi separation`= mean(near_separated, na.rm = TRUE),
    `MLE non-convergence`  = mean(mle_nonconv,    na.rm = TRUE),
    `Firth non-convergence`= mean(firth_nonconv,  na.rm = TRUE)
  ), by = x_var]
  setnames(rates, x_var, "x")

  rates_long <- melt(rates, id.vars = "x",
                     variable.name = "rate_type", value.name = "rate")
  rates_long[, rate_type := factor(rate_type,
    levels = c("Complete separation", "Near/quasi separation",
               "MLE non-convergence", "Firth non-convergence"))]

  p1 <- ggplot(rates_long, aes(x = x, y = rate, colour = rate_type, group = rate_type)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    x_scale +
    scale_colour_manual(values = c(
      "Complete separation"   = "#2c7bb6",
      "Near/quasi separation" = "#abd9e9",
      "MLE non-convergence"  = "#d7191c",
      "Firth non-convergence"= "#1a9641"
    )) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
    labs(x = x_label, y = "% of simulations", colour = NULL) +
    theme_bw(base_size = 13) +
    theme(legend.position = "top", panel.grid = element_blank()) +
    guides(colour = guide_legend(nrow = 2))

  # Plot 2: bias — MLE in non-separated runs only, Firth in all runs
  mle_nonsep <- dt[mle_nonconv == FALSE, .(x = get(x_var), estimate = mle_est,   method = "MLE")]
  firth_all  <- dt[firth_nonconv == FALSE, .(x = get(x_var), estimate = firth_est, method = "Firth")]
  est_long   <- rbindlist(list(mle_nonsep, firth_all))
  est_long[, method := factor(method, levels = c("MLE", "Firth"))]
  est_long[, bias   := estimate - true_beta]


  p2 <- ggplot(est_long, aes(x = x, y = bias, fill = method,
                              group = interaction(x, method))) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_boxplot() +
    x_scale +
    scale_fill_manual(values = c("MLE" = "#E07B54", "Firth" = "#5B8DB8")) +
    labs(
      x       = x_label,
      y       = expression("Bias = "~hat(beta) - beta),
      fill    = "Method"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "top", panel.grid = element_blank())

  # Plot 3: SE — MLE in non-separated runs only, Firth in all runs
  mle_se_dt   <- dt[mle_nonconv == FALSE, .(x = get(x_var), se = mle_se,   method = "MLE")]
  firth_se_dt <- dt[firth_nonconv == FALSE, .(x = get(x_var), se = firth_se, method = "Firth")]
  se_long     <- rbindlist(list(mle_se_dt, firth_se_dt))
  se_long[, method := factor(method, levels = c("MLE", "Firth"))]

  se_clip_hi <- quantile(se_long$se, 0.95, na.rm = TRUE)

  p3 <- ggplot(se_long, aes(x = x, y = se, fill = method,
                             group = interaction(x, method))) +
    geom_boxplot() +
    x_scale +
    scale_fill_manual(values = c("MLE" = "#E07B54", "Firth" = "#5B8DB8")) +
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
