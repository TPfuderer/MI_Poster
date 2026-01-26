# ============================================================
# Mini sim study (proof of concept) for miceRanger under MAR
# DGP: x1 -> x2 -> x3 (simple chained regression DGP)
# Missingness: MAR via logit model for Pr(x3 missing | x1,x2)
#
# Outputs:
# - RMSE on missing x3
# - Bias of mean(x3) over missing positions
# - Optional: quick regression coefficient bias (y ~ x1+x2+x3)
# - Optional: miceRanger diagnostic plots
# ============================================================

# ---- Packages ----
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  data.table,
  miceRanger,
  dplyr
)

set.seed(1)

# ============================================================
# 1) Data Generating Process (DGP)
#    (simple chained linear DGP; keep it small & transparent)
# ============================================================
simulate_complete_data <- function(n) {
  x1 <- rnorm(n, mean = 5, sd = 2)
  x2 <- 1 + 0.5 * x1 + rnorm(n, mean = 0, sd = 1)
  x3 <- 1 + 0.5 * x2 + rnorm(n, mean = 0, sd = 1)
  
  # Optional analysis outcome (for a quick "downstream" check)
  y <- 2 + 1.2 * x1 - 0.7 * x2 + 0.9 * x3 + rnorm(n, 0, 1)
  
  data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)
}

# ============================================================
# 2) MAR mechanism: probabilistic logit missingness for x3
#    Depends ONLY on observed covariates (x1, x2) -> MAR
# ============================================================
make_mar_missing_x3 <- function(dat, target_miss_rate = 0.30) {
  # Linear predictor (tweakable)
  # Higher x1 increases missingness; higher x2 decreases it
  eta <- -1 + 0.25 * dat$x1 - 0.35 * dat$x2
  
  # Calibrate intercept to hit approximate target missingness
  # (KISS calibration: shift eta by a constant found by 1D search)
  f <- function(shift) mean(stats::plogis(eta + shift)) - target_miss_rate
  shift_hat <- uniroot(f, interval = c(-10, 10))$root
  
  p_miss <- stats::plogis(eta + shift_hat)
  z_miss <- rbinom(nrow(dat), size = 1, prob = p_miss)
  
  dat_amp <- dat
  dat_amp$x3[z_miss == 1] <- NA
  
  list(
    data_incomplete = dat_amp,
    miss_indicator = z_miss,
    miss_prob = p_miss
  )
}

# ============================================================
# 3) Fit miceRanger and return completed datasets
# ============================================================
run_miceranger <- function(dat_incomplete, m = 5, maxiter = 3,
                           valueSelector = "value",
                           meanMatchCandidates = 0,
                           num.trees = 200,
                           returnModels = TRUE,
                           verbose = FALSE) {
  
  # miceRanger passes ... to ranger(), so we can set num.trees etc.
  miceObj <- miceRanger::miceRanger(
    dat_incomplete,
    m = m,
    maxiter = maxiter,
    valueSelector = valueSelector,
    meanMatchCandidates = meanMatchCandidates,
    returnModels = returnModels,
    verbose = verbose,
    num.trees = num.trees
  )
  
  imputed_list <- miceRanger::completeData(miceObj)
  list(miceObj = miceObj, imputed_list = imputed_list)
}

# ============================================================
# 4) Metrics (KISS)
#    - RMSE of x3 over originally missing positions
#    - Bias of mean(x3) over missing positions
#    - Optional: quick regression coefficient bias
# ============================================================
rmse <- function(truth, est) sqrt(mean((truth - est)^2))

evaluate_imputations <- function(dat_true, dat_incomplete, imputed_list) {
  miss_idx <- which(is.na(dat_incomplete$x3))
  if (length(miss_idx) == 0) stop("No missing values in x3 - nothing to evaluate.")
  
  # For each imputed dataset, evaluate x3 at missing positions
  per_ds <- lapply(imputed_list, function(d) {
    x3_hat <- d$x3[miss_idx]
    x3_true <- dat_true$x3[miss_idx]
    
    list(
      rmse_x3 = rmse(x3_true, x3_hat),
      bias_mean_x3 = mean(x3_hat) - mean(x3_true),
      
      # Downstream quick check: regression coefficients
      # (This is *not* a full MI pooling workflow; just a sanity check.)
      beta_x3 = coef(lm(y ~ x1 + x2 + x3, data = d))["x3"]
    )
  })
  
  # Average across m imputations (simple average; KISS)
  rmse_avg <- mean(sapply(per_ds, `[[`, "rmse_x3"))
  bias_avg <- mean(sapply(per_ds, `[[`, "bias_mean_x3"))
  beta_x3_avg <- mean(sapply(per_ds, `[[`, "beta_x3"))
  
  # Truth for the downstream coefficient
  beta_x3_true <- coef(lm(y ~ x1 + x2 + x3, data = dat_true))["x3"]
  
  data.frame(
    rmse_x3 = rmse_avg,
    bias_mean_x3 = bias_avg,
    beta_x3_hat = beta_x3_avg,
    beta_x3_true = beta_x3_true,
    beta_x3_bias = beta_x3_avg - beta_x3_true
  )
}

# ============================================================
# 5) Optional diagnostics (built-in miceRanger plots)
#    See miceRanger diagnostic vignette / README for these.
# ============================================================
plot_miceranger_diagnostics <- function(miceObj) {
  # These functions are documented in the miceRanger diagnostics vignette. :contentReference[oaicite:1]{index=1}
  print(plotDistributions(miceObj, vars = "allNumeric"))
  print(plotCorrelations(miceObj, vars = "allNumeric"))
  print(plotVarConvergence(miceObj, vars = "allNumeric"))
  print(plotModelError(miceObj, vars = "allNumeric"))
  print(plotVarImportance(miceObj, vars = "allNumeric"))
  print(plotImputationVariance(miceObj, vars = "allNumeric", monteCarloSimulations = 200))
}

# ============================================================
# 6) Simulation driver (mini study)
# ============================================================
mini_sim <- function(
    R = 30,
    n = 1000,
    target_miss_rate = 0.30,
    m = 5,
    maxiter = 3,
    num.trees = 200,
    make_plots_once = TRUE
) {
  results <- vector("list", R)
  
  for (r in seq_len(R)) {
    dat_true <- simulate_complete_data(n)
    mar <- make_mar_missing_x3(dat_true, target_miss_rate = target_miss_rate)
    dat_inc <- mar$data_incomplete
    
    fit <- run_miceranger(
      dat_incomplete = dat_inc,
      m = m,
      maxiter = maxiter,
      num.trees = num.trees,
      returnModels = TRUE,
      verbose = FALSE
    )
    
    # Optional: show diagnostics for the first replication
    if (make_plots_once && r == 1) {
      message("Plotting miceRanger diagnostics for replication 1 ...")
      plot_miceranger_diagnostics(fit$miceObj)
    }
    
    res_r <- evaluate_imputations(dat_true, dat_inc, fit$imputed_list)
    res_r$rep <- r
    res_r$miss_rate_realized <- mean(is.na(dat_inc$x3))
    
    results[[r]] <- res_r
    message(sprintf("Rep %d/%d done | miss=%.3f | RMSE=%.4f",
                    r, R, res_r$miss_rate_realized, res_r$rmse_x3))
  }
  
  dplyr::bind_rows(results)
}

# ============================================================
# 7) Run the mini sim
# ============================================================
out <- mini_sim(
  R = 20,              # keep small for proof-of-concept
  n = 1500,            # enough to stabilize a bit
  target_miss_rate = 0.30,
  m = 5,
  maxiter = 3,
  num.trees = 200,
  make_plots_once = TRUE
)

print(out)

summary_tab <- out %>%
  summarise(
    RMSE_x3_mean = mean(rmse_x3),
    RMSE_x3_sd   = sd(rmse_x3),
    BiasMean_x3_mean = mean(bias_mean_x3),
    Beta_x3_bias_mean = mean(beta_x3_bias),
    Beta_x3_bias_sd   = sd(beta_x3_bias),
    miss_rate_mean = mean(miss_rate_realized)
  )

print(summary_tab)

#8B Properness diagnostics (this is the key part)
set.seed(123)

dat_true <- simulate_complete_data(n = 1500)
mar <- make_mar_missing_x3(dat_true, target_miss_rate = 0.30)
dat_inc <- mar$data_incomplete

fit <- run_miceranger(
  dat_incomplete = dat_inc,
  m = 20,
  maxiter = 30,
  num.trees = 250,
  returnModels = TRUE
)


check_between_variance <- function(imputed_list) {
  means <- sapply(imputed_list, function(d) mean(d$x3, na.rm = TRUE))
  var(means)
}

check_between_variance(fit$imputed_list)

#
mc_error <- function(imputed_list) {
  means <- sapply(imputed_list, function(d) mean(d$x3, na.rm = TRUE))
  sd(means) / sqrt(length(means))
}

mc_error(fit$imputed_list)
#
ppc_x3_mean <- function(dat_true, imputed_list) {
  
  mu_true <- mean(dat_true$x3)
  mu_imp  <- sapply(imputed_list, function(d) mean(d$x3))
  
  hist(mu_imp,
       breaks = 10,
       col = "lightblue",
       main = "Posterior predictive check: mean(x3)",
       xlab = "Mean x3 (imputed)")
  abline(v = mu_true, col = "red", lwd = 2)
}

plotVarConvergence(fit$miceObj)
plotImputationVariance(fit$miceObj)
plotModelError(fit$miceObj)
plotVarImportance(fit$miceObj)
