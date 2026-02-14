# ============================================================
# BURGETTE & REITER STYLE SIMULATION (PROPERNESS)
# ============================================================
#### Load previous 
#sim_results <- readRDS(
  #"C:/Users/pfudi/PycharmProjects/MI_Poster/tristan stuff test/burgette_sim_2026-02-13.rds"
#)


#coef_mat <- sim_results$coef_mat
#se_mat   <- sim_results$se_mat
#Imp_store <- sim_results$Imp_store
##################

################## 
## Crash save 
#latest <- tail(sort(list.files(pattern="checkpoint_run_.*\\.rds$")),1)
#chk <- readRDS(latest)

#coef_store <- chk$coef_store
#se_store   <- chk$se_store
#start_run  <- chk$last_run + 1
###################

pacman::p_load(miceRanger, tidyverse, MASS, mice)

set.seed(123)

# -----------------------------
# Parameters 1:1 Burgette
# -----------------------------
n  <- 1000
p  <- 10
R  <- 15
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
rmse_store <- vector("list", R)
cor_store  <- numeric(R)


form_true <- Y ~ X1 + X2 + X3 + X8 + X9 +
  I(X3^2) + X1:X2 + X8:X9
options(ranger.num.threads = 4)

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
  # ---------------------------------
  # EXTRA ML DIAGNOSTICS
  # ---------------------------------
  ######## ADDD IIINNN LOOOOPPP
  # Use first imputed dataset for ML evaluation
  imp1 <- imputed_list[[1]]
  
  # Store true full data before missingness
  true_full <- full_df
  
  # Identify missing cells
  miss_mask <- is.na(miss_df)
  
  # 1️⃣ RMSE for missing entries only
  rmse_vals <- sapply(colnames(miss_df), function(v) {
    idx <- which(is.na(miss_df[[v]]))
    if(length(idx) == 0) return(NA_real_)
    sqrt(mean(
      (imp1[[v]][idx] - true_full[[v]][idx])^2
    ))
  })
  
  
  # 2️⃣ Correlation matrix recovery (Frobenius norm)
  cor_true <- cor(true_full)
  cor_imp  <- cor(imp1)
  
  cor_error <- sqrt(sum((cor_true - cor_imp)^2))
  
  # Store
  rmse_store[[r]] <- rmse_vals
  cor_store[r]    <- cor_error
  
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
  
  # ---- CHECKPOINT SAVE EVERY 5 RUNS ----
  if (r %% 5 == 0) {
    
    sim_checkpoint <- list(
      coef_store = coef_store,
      se_store   = se_store,
      last_run   = r
    )
    
    file_name <- paste0("checkpoint_run_", r, ".rds")
    saveRDS(sim_checkpoint, file_name)
    
    saveRDS(sim_checkpoint, "checkpoint_latest_backup.rds")
    
    files <- sort(list.files(pattern = "checkpoint_run_.*\\.rds$"))
    
    if (length(files) > 2) {
      file.remove(files[1])
    }
    
    cat("Checkpoint saved at run", r, "\n")
  }
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

par(mfrow = c(2,2))
for(i in 1:4){
  plotDistributions(Imp_store[[i]], vars = "Y")
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
##################################
# Convert RMSE list to matrix
rmse_mat <- do.call(rbind, rmse_store)

mean_rmse <- colMeans(rmse_mat, na.rm = TRUE)
mean_cor_error <- mean(cor_store)

mean_rmse
mean_cor_error

bias_df <- data.frame(
  Parameter = names(bias),
  Bias = as.numeric(bias)
)

library(ggplot2)

ggplot(bias_df, aes(x = reorder(Parameter, Bias), y = Bias)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Bias of MI Estimators",
       y = "Bias",
       x = "")

cov_df <- data.frame(
  Parameter = names(coverage),
  Coverage = as.numeric(coverage)
)

ggplot(cov_df, aes(x = reorder(Parameter, Coverage), y = Coverage)) +
  geom_col() +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_minimal() +
  labs(title = "95% CI Coverage",
       y = "Coverage",
       x = "")
var_df <- data.frame(
  Parameter = names(emp_var),
  Empirical = emp_var,
  Rubin = mean_rubin_var
)

var_long <- tidyr::pivot_longer(
  var_df,
  cols = c("Empirical", "Rubin"),
  names_to = "Type",
  values_to = "Variance"
)

ggplot(var_long, aes(x = Parameter, y = Variance, fill = Type)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Empirical vs Rubin Variance")
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
  geom_point(data = data.frame(
    Parameter = names(true_beta),
    True = true_beta
  ),
  aes(x = Parameter, y = True),
  color = "red",
  size = 3) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Distribution of Estimates Across Simulations")

# # ============================================================
# # FINAL SAVE – FULL REPRODUCIBLE SIMULATION STATE
# # (only save you need)
# # ============================================================
diag_runs <- seq_along(Imp_store)

final_results <- list(
  # results
  coef_mat        = coef_mat,
  se_mat          = se_mat,
  bias            = bias,
  emp_var         = emp_var,
  mean_rubin_var  = mean_rubin_var,
  coverage        = coverage,

  # ML diagnostics
  rmse_store      = rmse_store,
  cor_store       = cor_store,
  mean_rmse       = mean_rmse,
  mean_cor_error  = mean_cor_error,

  # stored imputation objects for diagnostics
  diag_runs       = diag_runs,
  Imp_store       = Imp_store,

  # DGP/meta
  beta            = beta,
  true_beta       = true_beta,
  Sigma           = Sigma,
  form_true       = deparse(form_true),
  n               = n,
  p               = p,
  R               = R,
  m_val           = m_val,
  missing_mechanism = "MAR: logit(0.5*X9 - 0.5*X10)",

  # reproducibility
  seed            = 123,
  RNGkind         = RNGkind(),
  sessionInfo     = sessionInfo(),
  date            = Sys.time()
)

saveRDS(final_results, paste0("burgette_full_sim_", Sys.Date(), ".rds"))
cat("Saved: burgette_full_sim_", Sys.Date(), ".rds\n", sep = "")



