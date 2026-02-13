# ============================================================
# MICERANGER SIMULATION STUDY
# ============================================================

#ADDDDD MISSINGNESSS PROCESSSSSS


pacman::p_load(miceRanger, car, tidyverse)

set.seed(123)

# Full data (no missing)
full_df <- read_csv("tristan stuff test/FULL_real_dataset.csv")

R <- 200        # number of simulation repetitions
m_val <- 25     # imputations

# storage
coef_store <- list()
se_store   <- list()

for (r in 1:R) {
  
  cat("Simulation:", r, "\n")
  
  # ------------------------------------------------------------
  # 1) Impose missingness (example: 50% MCAR on X1–X8 + Y)
  # ------------------------------------------------------------
  
  miss_df <- full_df
  
  for (var in c("Y", paste0("X", 1:8))) {
    miss_index <- sample(1:nrow(miss_df), 
                         size = floor(0.5 * nrow(miss_df)))
    miss_df[miss_index, var] <- NA
  }
  
  # ------------------------------------------------------------
  # 2) Imputation
  # ------------------------------------------------------------
  
  Imp_RF <- miceRanger(
    miss_df,
    m = m_val,
    maxit = 10,
    num.threads = 4,
    verbose = FALSE
  )
  
  imputed_list <- completeData(Imp_RF)
  
  # ------------------------------------------------------------
  # 3) Model fitting
  # ------------------------------------------------------------
  
  yfit_list <- lapply(
    imputed_list,
    function(dat) lm(Y ~ X1 + X2 + X3 + X4 + X5 + 
                       X6 + X7 + X8 + X9 + X10,
                     data = dat)
  )
  
  class(yfit_list) <- "mira"
  pooled_rf <- mice::pool(yfit_list)
  
  pooled_table <- pooled_rf$pooled
  
  # store
  coef_store[[r]] <- pooled_table$estimate
  se_store[[r]]   <- sqrt(pooled_table$t)
}

# ============================================================
# POST-SIMULATION EVALUATION
# ============================================================

coef_mat <- do.call(rbind, coef_store)
se_mat   <- do.call(rbind, se_store)

# Bias (if you know true beta)
true_beta <- coef(lm(Y ~ X1 + X2 + X3 + X4 + X5 + 
                       X6 + X7 + X8 + X9 + X10,
                     data = full_df))

bias <- colMeans(coef_mat) - true_beta

# Empirical variance
emp_var <- apply(coef_mat, 2, var)

# Mean Rubin variance
mean_rubin_var <- colMeans(se_mat^2)

bias
emp_var
mean_rubin_var
