# ============================================================
# TEST A1 (Standalone): 0% missing baseline inference
# Goal: check if analysis model + SE are calibrated without MI
# ============================================================

pacman::p_load(tidyverse, MASS)

set.seed(123)

# -----------------------------
# Parameters (FAST DEBUG)
# -----------------------------
n <- 1000
p <- 10
R <- 10

beta <- c(
  b0 = 0,
  b1 = 0.5,
  b2 = 0.5,
  b3 = 0.5,
  b4 = 0.5,
  b5 = 0.5,
  b6 = 1,
  b7 = 1,
  b8 = 1
)

Sigma <- diag(p)
Sigma[1:4, 1:4]   <- 0.5; diag(Sigma[1:4, 1:4])   <- 1
Sigma[5:10, 5:10] <- 0.3; diag(Sigma[5:10, 5:10]) <- 1

form_true <- Y ~ X1 + X2 + X3 + X8 + X9 + I(X3^2) + X1:X2 + X8:X9

true_beta <- c(
  "(Intercept)" = unname(beta["b0"]),
  "X1"          = unname(beta["b1"]),
  "X2"          = unname(beta["b2"]),
  "X3"          = unname(beta["b3"]),
  "X8"          = unname(beta["b4"]),
  "X9"          = unname(beta["b5"]),
  "I(X3^2)"     = unname(beta["b6"]),
  "X1:X2"       = unname(beta["b7"]),
  "X8:X9"       = unname(beta["b8"])
)


# -----------------------------
# Storage
# -----------------------------
theta_hat <- matrix(NA_real_, nrow = R, ncol = length(true_beta))
se_hat    <- matrix(NA_real_, nrow = R, ncol = length(true_beta))
colnames(theta_hat) <- names(true_beta)
colnames(se_hat)    <- names(true_beta)

# ============================================================
# SIM LOOP (NO MISSINGNESS)
# ============================================================

for (r in 1:R) {
  
  # 1) Generate complete DGP
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
  
  # 2) Fit analysis model directly
  fit <- lm(form_true, data = full_df)
  coefs <- coef(fit)
  ses   <- summary(fit)$coef[, "Std. Error"]
  
  # 3) Store only the coefficients we care about
  theta_hat[r, ] <- coefs[names(true_beta)]
  se_hat[r, ]    <- ses[names(true_beta)]
}

# ============================================================
# OUTPUT: bias_full, cover_full, sd(z_full)
# ============================================================

bias_full <- colMeans(theta_hat) - true_beta

cover_full <- colMeans(
  (theta_hat - 1.96 * se_hat <= rep(true_beta, each = R)) &
    (theta_hat + 1.96 * se_hat >= rep(true_beta, each = R))
)

z_full <- sweep(theta_hat, 2, true_beta, "-") / se_hat
sd_z_full <- apply(z_full, 2, sd)
mean_z_full <- apply(z_full, 2, mean)

out <- data.frame(
  param       = names(true_beta),
  true_beta   = as.numeric(true_beta),
  bias_full   = as.numeric(bias_full),
  cover_full  = as.numeric(cover_full),
  mean_z_full = as.numeric(mean_z_full),
  sd_z_full   = as.numeric(sd_z_full)
) %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))

cat("\n=====================================================\n")
cat("TEST A1: 0% missing baseline inference\n")
cat("R =", R, "| n =", n, "\n")
cat("Targets: bias ~ 0, cover ~ 0.95, mean(z) ~ 0, sd(z) ~ 1\n")
cat("=====================================================\n\n")

print(out)

cat("\nQuick summary:\n")
cat("Mean abs(bias):", round(mean(abs(bias_full)), 4), "\n")
cat("Mean coverage:", round(mean(cover_full), 4), "\n")
cat("Mean sd(z):", round(mean(sd_z_full), 4), "\n")

