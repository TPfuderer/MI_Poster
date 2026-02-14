# ============================================================
# BURGETTE & REITER STYLE SIMULATION (PROPERNESS) — CLEAN FILE
# ============================================================

pacman::p_load(miceRanger, tidyverse, MASS, mice, ggplot2)

set.seed(123)

# -----------------------------
# Parameters 1:1 Burgette
# -----------------------------
n     <- 1000
p     <- 10
R     <- 100         # set to 15 for quick tests
m_val <- 20

beta <- c(
  b0 = 0,
  b1 = 0.5,
  b2 = 0.5,
  b3 = 0.5,
  b4 = 0.5,
  b5 = 0.5,
  b6 = 1,    # X3^2
  b7 = 1,    # X1*X2
  b8 = 1     # X8*X9
)

Sigma <- diag(p)
Sigma[1:4, 1:4]     <- 0.5; diag(Sigma[1:4, 1:4])     <- 1
Sigma[5:10, 5:10]   <- 0.3; diag(Sigma[5:10, 5:10])   <- 1

logit_p <- function(z) 1 / (1 + exp(-z))

form_true <- Y ~ X1 + X2 + X3 + X8 + X9 + I(X3^2) + X1:X2 + X8:X9

true_beta <- c(
  "(Intercept)" = beta["b0"],
  "X1"          = beta["b1"],
  "X2"          = beta["b2"],
  "X3"          = beta["b3"],
  "X8"          = beta["b4"],
  "X9"          = beta["b5"],
  "I(X3^2)"     = beta["b6"],
  "X1:X2"       = beta["b7"],
  "X8:X9"       = beta["b8"]
)

options(ranger.num.threads = 4)

# strg + alt + c => to uncomment!
# ============================================================
# (A) SIMULATION LOOP — TEST VERSION (R = 5)
# ============================================================

# Storage
coef_store <- vector("list", R)
se_store   <- vector("list", R)
rmse_store <- vector("list", R)
cor_store  <- numeric(R)

# For test run: just store ALL runs (since R=5)
diag_runs <- 1:R
Imp_store <- vector("list", R)

for (r in 1:R) {
  
  cat("Simulation:", r, "\n")
  
  # ------------------------------------------------------------
  # 1) Generate NEW DGP
  # ------------------------------------------------------------
  X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(X) <- paste0("X", 1:10)
  
  epsilon <- rnorm(n)
  
  Y <- beta["b0"] +
    beta["b1"] * X[, "X1"] +
    beta["b2"] * X[, "X2"] +
    beta["b3"] * X[, "X3"] +
    beta["b4"] * X[, "X8"] +
    beta["b5"] * X[, "X9"] +
    beta["b6"] * X[, "X3"]^2 +
    beta["b7"] * X[, "X1"] * X[, "X2"] +
    beta["b8"] * X[, "X8"] * X[, "X9"] +
    epsilon
  
  full_df <- data.frame(Y, X)
  
  # ------------------------------------------------------------
  # 2) MAR Missingness
  # ------------------------------------------------------------
  miss_df <- full_df
  p_miss  <- logit_p(0.5 * miss_df$X9 - 0.5 * miss_df$X10)
  
  for (v in c("Y", paste0("X", 1:8))) {
    miss_indicator <- rbinom(n, 1, p_miss)
    miss_df[miss_indicator == 1, v] <- NA
  }
  
  # ------------------------------------------------------------
  # 3) Imputation
  # ------------------------------------------------------------
  Imp_RF <- miceRanger(
    miss_df,
    m = m_val,
    maxit = 10,
    num.threads = 4,
    verbose = FALSE
  )
  
  # Store ALL runs (since R=5)
  Imp_store[[r]] <- Imp_RF
  
  imputed_list <- completeData(Imp_RF)
  
  # ------------------------------------------------------------
  # 3b) ML Diagnostics
  # ------------------------------------------------------------
  imp1      <- as.data.frame(imputed_list[[1]])
  true_full <- full_df
  
  rmse_vals <- sapply(colnames(miss_df), function(v) {
    idx <- which(is.na(miss_df[[v]]))
    if (length(idx) == 0) return(NA_real_)
    sqrt(mean((imp1[[v]][idx] - true_full[[v]][idx])^2))
  })
  
  cor_true  <- cor(true_full)
  cor_imp   <- cor(imp1)
  cor_error <- sqrt(sum((cor_true - cor_imp)^2))
  
  rmse_store[[r]] <- rmse_vals
  cor_store[r]    <- cor_error
  
  # ------------------------------------------------------------
  # 4) Nonlinear Analysis Model
  # ------------------------------------------------------------
  yfit_list <- lapply(imputed_list, function(dat)
    lm(form_true, data = dat)
  )
  
  class(yfit_list) <- "mira"
  pooled_rf <- mice::pool(yfit_list)
  pooled_table <- pooled_rf$pooled
  
  coef_store[[r]] <- pooled_table$estimate
  se_store[[r]]   <- sqrt(pooled_table$t)
}

# ============================================================
# EVALUATION
# ============================================================

coef_mat <- do.call(rbind, coef_store)
se_mat   <- do.call(rbind, se_store)

colnames(coef_mat) <- names(true_beta)
colnames(se_mat)   <- names(true_beta)

bias <- colMeans(coef_mat) - true_beta
emp_var <- apply(coef_mat, 2, var)
mean_rubin_var <- colMeans(se_mat^2)

coverage <- colMeans(
  (coef_mat - 1.96 * se_mat <= true_beta) &
    (coef_mat + 1.96 * se_mat >= true_beta)
)

rmse_mat <- do.call(rbind, rmse_store)
mean_rmse <- colMeans(rmse_mat, na.rm = TRUE)
mean_cor_error <- mean(cor_store)

# ============================================================
# FINAL SAVE
# ============================================================

final_results <- list(
  coef_mat        = coef_mat,
  se_mat          = se_mat,
  bias            = bias,
  emp_var         = emp_var,
  mean_rubin_var  = mean_rubin_var,
  coverage        = coverage,
  rmse_store      = rmse_store,
  cor_store       = cor_store,
  mean_rmse       = mean_rmse,
  mean_cor_error  = mean_cor_error,
  diag_runs       = diag_runs,
  Imp_store       = Imp_store,
  beta            = beta,
  true_beta       = true_beta,
  Sigma           = Sigma,
  form_true       = deparse(form_true),
  n               = n,
  p               = p,
  R               = R,
  m_val           = m_val,
  missing_mechanism = "MAR: logit(0.5*X9 - 0.5*X10)",
  seed            = 123,
  RNGkind         = RNGkind(),
  sessionInfo     = sessionInfo(),
  date            = Sys.time()
)

saveRDS(final_results, "burgette_full_sim_TEST.rds")

cat("\nTest simulation complete and saved.\n")

# ============================================================
# (B) DEFAULT BEHAVIOR: LOAD SAVED RESULTS + RUN DIAGNOSTICS
# ============================================================

# --- set your file path here ---
# Example:
# sim_path <- "C:/Users/pfudi/PycharmProjects/MI_Poster/tristan stuff test/burgette_full_sim_2026-02-14.rds"
sim_path <- "burgette_full_sim_2026-02-14.rds"

sim <- readRDS(sim_path)
 
coef_mat       <- sim$coef_mat
se_mat         <- sim$se_mat
bias           <- sim$bias
emp_var        <- sim$emp_var
mean_rubin_var <- sim$mean_rubin_var
coverage       <- sim$coverage

rmse_store     <- sim$rmse_store
cor_store      <- sim$cor_store
mean_rmse      <- sim$mean_rmse
mean_cor_error <- sim$mean_cor_error

Imp_store      <- sim$Imp_store
diag_runs      <- sim$diag_runs

true_beta      <- sim$true_beta

# ============================================================
# (C1) NORMAL IMPUTATION MODEL DIAGNOSTICS (4 stored runs)
# ============================================================

library(gridExtra)

vars_diag <- "Y"

cat("\nStored diagnostic runs:", paste(diag_runs, collapse = ", "), "\n\n")

# ------------------------------------------------------------
# Helper to safely extract miceRanger plot grob
# ------------------------------------------------------------

extract_plot <- function(plot_obj) {
  if (is.list(plot_obj) && length(plot_obj) == 1) {
    return(plot_obj[[1]])
  }
  return(plot_obj)
}

# ------------------------------------------------------------
# 1️⃣ Distributions
# ------------------------------------------------------------

dist_plots <- lapply(Imp_store[1:4], function(obj) {
  extract_plot(plotDistributions(obj, vars = vars_diag))
})

do.call(grid.arrange, c(dist_plots, ncol = 2))

# ------------------------------------------------------------
# 2️⃣ Convergence
# ------------------------------------------------------------

conv_plots <- lapply(Imp_store[1:4], function(obj) {
  extract_plot(plotVarConvergence(obj, vars = vars_diag))
})

do.call(grid.arrange, c(conv_plots, ncol = 2))

# ------------------------------------------------------------
# 3️⃣ Imputation Variance
# ------------------------------------------------------------

var_plots <- lapply(Imp_store[1:4], function(obj) {
  extract_plot(plotImputationVariance(obj, vars = vars_diag))
})

do.call(grid.arrange, c(var_plots, ncol = 2))

# ------------------------------------------------------------
# 4️⃣ RF Model Error
# ------------------------------------------------------------

error_plots <- lapply(Imp_store[1:4], function(obj) {
  extract_plot(plotModelError(obj))
})

do.call(grid.arrange, c(error_plots, ncol = 2))

# ------------------------------------------------------------
# ML Quality Summary
# ------------------------------------------------------------

cat("\nMean RMSE (missing entries only):\n")
print(mean_rmse)

cat("\nMean correlation-matrix error (Frobenius norm):\n")
print(mean_cor_error)


# ============================================================
# (C2) PROPERNESS DIAGNOSTICS (bias/coverage/variance/estimates)
# ============================================================

# --- Bias plot ---
bias_df <- data.frame(Parameter = names(bias), Bias = as.numeric(bias))

ggplot(bias_df, aes(x = reorder(Parameter, Bias), y = Bias)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Bias of MI Estimators", x = "", y = "Bias")

# --- Coverage plot ---
cov_df <- data.frame(Parameter = names(coverage), Coverage = as.numeric(coverage))

ggplot(cov_df, aes(x = reorder(Parameter, Coverage), y = Coverage)) +
  geom_col() +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_minimal() +
  labs(title = "95% CI Coverage", x = "", y = "Coverage")

# --- Empirical vs Rubin variance ---
var_df <- data.frame(
  Parameter = names(emp_var),
  Empirical = as.numeric(emp_var),
  Rubin     = as.numeric(mean_rubin_var)
)

var_long <- tidyr::pivot_longer(
  var_df,
  cols = c("Empirical", "Rubin"),
  names_to = "Type",
  values_to = "Variance"
)

ggplot(var_long, aes(x = reorder(Parameter, Variance), y = Variance, fill = Type)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Empirical vs Rubin Variance", x = "", y = "Variance")

# --- Distribution of estimates across simulations (boxplot + true values) ---
coef_df <- as.data.frame(coef_mat)
coef_df$Simulation <- 1:nrow(coef_df)

coef_long <- tidyr::pivot_longer(
  coef_df,
  cols = -Simulation,
  names_to = "Parameter",
  values_to = "Estimate"
)

ggplot(coef_long, aes(x = Parameter, y = Estimate)) +
  geom_boxplot() +
  geom_point(
    data = data.frame(Parameter = names(true_beta), True = as.numeric(true_beta)),
    aes(x = Parameter, y = True),
    color = "red",
    size = 3
  ) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Distribution of Estimates Across Simulations")

# ------------------------------------------------------------
# 0) If loading from saved results:
# ------------------------------------------------------------
# final_results <- readRDS("burgette_full_sim_YYYY-MM-DD.rds")
# attach values if needed:
# coef_mat       <- final_results$coef_mat
# se_mat         <- final_results$se_mat
# bias           <- final_results$bias
# emp_var        <- final_results$emp_var
# mean_rubin_var <- final_results$mean_rubin_var
# coverage       <- final_results$coverage
# true_beta      <- final_results$true_beta

library(ggplot2)
library(tidyr)
library(dplyr)

# ------------------------------------------------------------
# 1️⃣ Bias Distribution Across Simulations (SAE style)
# ------------------------------------------------------------

coef_df <- as.data.frame(coef_mat)
coef_df$Simulation <- 1:nrow(coef_df)

coef_long <- pivot_longer(
  coef_df,
  cols = -Simulation,
  names_to = "Parameter",
  values_to = "Estimate"
)

true_df <- data.frame(
  Parameter = names(true_beta),
  True = as.numeric(true_beta)
)

ggplot(coef_long, aes(x = Parameter, y = Estimate)) +
  geom_boxplot(fill = "lightblue") +
  geom_point(data = true_df,
             aes(x = Parameter, y = True),
             color = "red",
             size = 3) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Distribution of Estimates (SAE Simulation Style)",
       subtitle = "Red dot = true parameter")

# ------------------------------------------------------------
# 2️⃣ Bias Bar Plot
# ------------------------------------------------------------

bias_df <- data.frame(
  Parameter = names(bias),
  Bias = as.numeric(bias)
)

ggplot(bias_df, aes(x = reorder(Parameter, Bias), y = Bias)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Bias of MI Estimators",
       y = "Bias")

# ------------------------------------------------------------
# 3️⃣ Coverage Plot (Properness Diagnostic)
# ------------------------------------------------------------

cov_df <- data.frame(
  Parameter = names(coverage),
  Coverage = as.numeric(coverage)
)

ggplot(cov_df, aes(x = reorder(Parameter, Coverage), y = Coverage)) +
  geom_col(fill = "darkorange") +
  geom_hline(yintercept = 0.95,
             linetype = "dashed",
             color = "red") +
  coord_flip() +
  theme_minimal() +
  labs(title = "95% CI Coverage",
       subtitle = "Dashed line = nominal 0.95",
       y = "Coverage")

# ------------------------------------------------------------
# 4️⃣ Variance Calibration Plot
#    (Empirical vs Rubin Variance)
# ------------------------------------------------------------

var_df <- data.frame(
  Parameter = names(emp_var),
  Empirical = emp_var,
  Rubin = mean_rubin_var
)

var_long <- pivot_longer(
  var_df,
  cols = c("Empirical", "Rubin"),
  names_to = "Type",
  values_to = "Variance"
)

ggplot(var_long, aes(x = Parameter, y = Variance, fill = Type)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Empirical vs Rubin Variance",
       subtitle = "Proper MI ⇒ Rubin ≈ Empirical")

# ------------------------------------------------------------
# 5️⃣ Interval Length Diagnostic
# ------------------------------------------------------------

ci_length <- 2 * 1.96 * se_mat
mean_ci_length <- colMeans(ci_length)

ci_df <- data.frame(
  Parameter = names(mean_ci_length),
  Mean_CI_Length = mean_ci_length
)

ggplot(ci_df, aes(x = reorder(Parameter, Mean_CI_Length),
                  y = Mean_CI_Length)) +
  geom_col(fill = "purple") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Average 95% CI Length",
       subtitle = "Detects over/under confidence")

# EXTRAS:
library(tidyverse)

plot_missing_vs_driver <- function(miss_df, driver = "X9", v = "Y", bins = 20) {
  df <- miss_df %>%
    transmute(
      driver = .data[[driver]],
      miss   = as.integer(is.na(.data[[v]]))
    ) %>%
    mutate(bin = ntile(driver, bins)) %>%
    group_by(bin) %>%
    summarise(
      driver_mean = mean(driver, na.rm = TRUE),
      miss_rate   = mean(miss),
      n           = n(),
      .groups = "drop"
    )
  
  ggplot(df, aes(driver_mean, miss_rate)) +
    geom_point() +
    geom_line() +
    labs(
      title = paste0("Missing rate of ", v, " vs ", driver),
      x = paste0(driver, " (binned mean)"),
      y = paste0("P(", v, " is missing)")
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1))
}

# Examples:
plot_missing_vs_driver(miss_df, driver="X9",  v="Y")
plot_missing_vs_driver(miss_df, driver="X10", v="Y")

