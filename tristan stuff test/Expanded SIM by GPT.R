# ============================================================
# BURGETTE & REITER STYLE SIMULATION (PROPERNESS)
# ============================================================

pacman::p_load(miceRanger, tidyverse, MASS, mice)

set.seed(123)

# -----------------------------
# Parameters 1:1 Burgette
# -----------------------------
n  <- 1000
p  <- 10
R  <- 100
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
Sigma[1:4,1:4]  <- 0.5; diag(Sigma[1:4,1:4])  <- 1
Sigma[5:10,5:10] <- 0.3; diag(Sigma[5:10,5:10]) <- 1

logit_p <- function(z) 1/(1+exp(-z))

coef_store <- list()
se_store   <- list()
Imp_store  <- vector("list",4)

form_true <- Y ~ X1 + X2 + X3 + X8 + X9 +
  I(X3^2) + X1:X2 + X8:X9

for(r in 1:R){
  
  cat("Simulation:", r, "\n")
  
  # ------------------------------------------------------------
  # 1) Generate NEW DGP each run
  # ------------------------------------------------------------
  X <- mvrnorm(n, mu = rep(0,p), Sigma = Sigma)
  colnames(X) <- paste0("X",1:10)
  
  epsilon <- rnorm(n)
  
  Y <- beta["b0"] +
    beta["b1"]*X[,1] +
    beta["b2"]*X[,2] +
    beta["b3"]*X[,3] +
    beta["b4"]*X[,8] +
    beta["b5"]*X[,9] +
    beta["b6"]*X[,3]^2 +
    beta["b7"]*X[,1]*X[,2] +
    beta["b8"]*X[,8]*X[,9] +
    epsilon
  
  full_df <- data.frame(Y, X)
  
  # ------------------------------------------------------------
  # 2) Impose MAR missingness
  # ------------------------------------------------------------
  miss_df <- full_df
  p_miss <- logit_p(0.5*miss_df$X9 - 0.5*miss_df$X10)
  
  for(v in c("Y", paste0("X",1:8))){
    miss_indicator <- rbinom(n,1,p_miss)
    miss_df[miss_indicator==1, v] <- NA
  }
  
  # ------------------------------------------------------------
  # 3) Imputation
  # ------------------------------------------------------------
  Imp_RF <- miceRanger(
    miss_df,
    m = m_val,
    maxit = 5,
    num.threads = 4,
    verbose = FALSE
  )
  
  if(r <= 4) Imp_store[[r]] <- Imp_RF
  
  imputed_list <- completeData(Imp_RF)
  
  # ------------------------------------------------------------
  # 4) Correct nonlinear analysis model
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

# True parameters from DGP
true_beta <- c(
  "(Intercept)" = beta["b0"],
  "X1" = beta["b1"],
  "X2" = beta["b2"],
  "X3" = beta["b3"],
  "X8" = beta["b4"],
  "X9" = beta["b5"],
  "I(X3^2)" = beta["b6"],
  "X1:X2" = beta["b7"],
  "X8:X9" = beta["b8"]
)

colnames(coef_mat) <- names(true_beta)
colnames(se_mat)   <- names(true_beta)

bias <- colMeans(coef_mat) - true_beta
emp_var <- apply(coef_mat,2,var)
mean_rubin_var <- colMeans(se_mat^2)

coverage <- colMeans(
  (coef_mat - 1.96*se_mat <= true_beta) &
    (coef_mat + 1.96*se_mat >= true_beta)
)

bias
emp_var
mean_rubin_var
coverage


# ============================================================
# DIAGNOSTICS: FIRST 4 RUNS (2x2)
# ============================================================

vars_diag <- "Y"

par(mfrow=c(2,2))
for(i in 1:4){
  plotDistributions(Imp_store[[i]], vars = vars_diag)
  mtext(paste("Run",i),3,1)
}

par(mfrow=c(2,2))
for(i in 1:4){
  plotVarConvergence(Imp_store[[i]], vars = vars_diag)
  mtext(paste("Run",i),3,1)
}

# ---- (C) Imputation variance / uncertainty (Run 1–4) ----
par(mfrow = c(2, 2), mar = c(3, 3, 3, 1))
for (i in 1:4) {
  plotImputationVariance(Imp_store[[i]], vars = vars_diag)
  mtext(paste0("Run ", i, " — plotImputationVariance"), side = 3, line = 1, cex = 0.9)
}

# ---- (D) RF model error (Run 1–4) ----
par(mfrow = c(2, 2), mar = c(3, 3, 3, 1))
for (i in 1:4) {
  plotModelError(Imp_store[[i]])
  mtext(paste0("Run ", i, " — plotModelError"), side = 3, line = 1, cex = 0.9)
}

