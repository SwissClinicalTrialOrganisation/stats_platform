library(brglm2)
library(logistf)
library(data.table)
library(parallel)

# Parallel streams of random numbers (reproducible for a fixed seed AND fixed
# number of cores; mclapply is fork-based, so it works on Linux/macOS/WSL only).
RNGkind("L'Ecuyer-CMRG")
set.seed(42)
n_cores <- as.integer(Sys.getenv("N_CORES", max(1, detectCores() - 1)))

true_beta          <- log(1/3)  # OR = 3
n_sims             <- 2000
conv_threshold     <- 10      # |estimate| or SE beyond this flags non-convergence
near_sep_min_count <- 5       # any 2x2 cell with fewer than this many subjects flags near-separation

# ── Helpers ───────────────────────────────────────────────────────────────────

intercept_for_prevalence <- function(target, true_beta) {
  uniroot(
    function(a) 0.5 * plogis(a) + 0.5 * plogis(a + true_beta) - target,
    interval = c(-20, 20)
  )$root
}

# Number of levels of x (0/1) in which the outcome is constant, i.e. the number
# of x groups that have an empty cell in the 2x2 table. Assumes y has both values.
#   2 -> complete separation:       x perfectly predicts y
#   1 -> quasi-complete separation: one empty cell, the other group has both outcomes
#   0 -> no separation
n_constant_levels <- function(x, y) {
  sum(sapply(unique(x), function(xi) length(unique(y[x == xi])) == 1))
}

# Near separation: at least one cell of the 2x2 table holds fewer than
# near_sep_min_count subjects (includes empty cells, i.e. also quasi-complete separation).
is_near_separated <- function(x, y, near_sep_min_count = 5) {
  cell_counts <- sapply(c(0, 1), function(xi)
    sapply(c(0, 1), function(yi) sum(x == xi & y == yi))
  )
  any(cell_counts < near_sep_min_count)
}

extract_coef_se <- function(fit) {
  if (is.null(fit)) return(c(est = NA_real_, se = NA_real_))
  s <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  if (is.null(s) || !("x" %in% rownames(s)))
    return(c(est = NA_real_, se = NA_real_))
  c(est = s["x", "Estimate"], se = s["x", "Std. Error"])
}

# logistf: coefficient and SE from the penalised likelihood fit (sqrt of diag of vcov)
extract_coef_se_logistf <- function(fit) {
  if (is.null(fit)) return(c(est = NA_real_, se = NA_real_))
  est <- tryCatch(unname(coef(fit)["x"]),                 error = function(e) NA_real_)
  se  <- tryCatch(unname(sqrt(diag(vcov(fit)))["x"]),     error = function(e) NA_real_)
  c(est = est, se = se)
}

is_nonconv <- function(est, se, threshold) {
  is.na(est) | is.na(se) | abs(est) > threshold | se > threshold
}

run_sim <- function(n, intercept, n_sims, true_beta, conv_threshold, near_sep_min_count) {
  results <- mclapply(seq_len(n_sims), function(i) {
    x <- rbinom(n, 1, 0.5)
    y <- rbinom(n, 1, plogis(intercept + true_beta * x))

    if (length(unique(y)) < 2) return(NULL)

    n_const        <- n_constant_levels(x, y)
    complete       <- n_const == 2
    quasi_complete <- n_const == 1
    near_separated <- is_near_separated(x, y, near_sep_min_count)

    fit_mle   <- tryCatch(glm(y ~ x, family = binomial),                      error = function(e) NULL)
    fit_firth <- tryCatch(glm(y ~ x, family = binomial, method = "brglmFit"), error = function(e) NULL)
    fit_lgf   <- tryCatch(suppressWarnings(logistf(y ~ x, data = data.frame(x = x, y = y))),
                          error = function(e) NULL)

    mle_cs   <- extract_coef_se(fit_mle)
    firth_cs <- extract_coef_se(fit_firth)
    lgf_cs   <- extract_coef_se_logistf(fit_lgf)

    data.frame(
      sim            = i,
      mle_est        = mle_cs[["est"]],
      mle_se         = mle_cs[["se"]],
      firth_est      = firth_cs[["est"]],
      firth_se       = firth_cs[["se"]],
      complete       = complete,
      quasi_complete = quasi_complete,
      near_separated = near_separated,
      mle_nonconv    = is_nonconv(mle_cs[["est"]],   mle_cs[["se"]],   conv_threshold),
      firth_nonconv  = is_nonconv(firth_cs[["est"]], firth_cs[["se"]], conv_threshold),
      logistf_est    = lgf_cs[["est"]],
      logistf_se     = lgf_cs[["se"]],
      logistf_nonconv = is_nonconv(lgf_cs[["est"]], lgf_cs[["se"]], conv_threshold)
    )
  }, mc.cores = n_cores, mc.set.seed = TRUE)
  rbindlist(results)
}

# ── Simulation 1: varying sample size, fixed prevalence (~5%) ────────────────
n_sizes    <- c(10, 20, 30, 50, 100, 200, 300, 500)
intercept1 <- intercept_for_prevalence(0.05, true_beta)

sim1 <- rbindlist(lapply(n_sizes, function(n) {
  cat(sprintf("[Sim 1] n = %d\n", n))
  res <- run_sim(n, intercept1, n_sims, true_beta, conv_threshold, near_sep_min_count)
  res$n <- n
  res
}))

# ── Simulation 2: fixed sample size (n=100), varying prevalence ──────────────
n_fixed     <- 100
prevalences <- c(0.01, 0.02, 0.03, 0.05, 0.10, 0.20, 0.30, 0.50)
intercepts2 <- sapply(prevalences, intercept_for_prevalence, true_beta = true_beta)

sim2 <- rbindlist(lapply(seq_along(prevalences), function(j) {
  cat(sprintf("[Sim 2] prevalence = %.0f%%\n", prevalences[j] * 100))
  res <- run_sim(n_fixed, intercepts2[j], n_sims, true_beta, conv_threshold, near_sep_min_count)
  res$prevalence <- prevalences[j] * 100
  res
}))

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(sim1, "firth/sim1.rds")
saveRDS(sim2, "firth/sim2.rds")
saveRDS(
  list(
    true_beta          = true_beta,
    n_sims             = n_sims,
    conv_threshold     = conv_threshold,
    near_sep_min_count = near_sep_min_count,
    n_sizes            = n_sizes,
    n_fixed            = n_fixed,
    prevalences        = prevalences
  ),
  "firth/sim_params.rds"
)
cat("Simulation data saved to firth/sim1.rds, firth/sim2.rds, and firth/sim_params.rds\n")
