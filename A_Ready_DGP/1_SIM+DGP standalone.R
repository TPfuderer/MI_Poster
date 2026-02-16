# ============================================================
# BURGETTE STYLE SIMULATION — MULTI METHOD VERSION
# ============================================================

pacman::p_load(miceRanger, tidyverse, MASS, mice, parallel)

set.seed(123)

missing_mech <- "MAR"   # "MAR" or "MCAR"
#pi_mcar <- 0.20         # target missing prob for MCAR
# -----------------------------
# Parameters
# -----------------------------
n     <- 1000
p     <- 10
R     <- 20
m_val <- 10

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

logit_p <- function(z) 1 / (1 + exp(-z))

form_true <- Y ~ X1 + X2 + X3 + X8 + X9 + X3_sq + X1_X2 + X8_X9

true_beta <- c(
  "(Intercept)" = unname(beta["b0"]),
  "X1"          = unname(beta["b1"]),
  "X2"          = unname(beta["b2"]),
  "X3"          = unname(beta["b3"]),
  "X8"          = unname(beta["b4"]),
  "X9"          = unname(beta["b5"]),
  "X3_sq"       = unname(beta["b6"]),
  "X1_X2"       = unname(beta["b7"]),
  "X8_X9"       = unname(beta["b8"])
)

options(ranger.num.threads = 5)

# ============================================================
# METHODS + STORAGE
# ============================================================

methods <- c("rf_ranger", "pmm", "cart") #, "rf_mice"

results <- lapply(methods, function(x) {
  list(
    coef_store = vector("list", R),
    se_store   = vector("list", R),
    df_store   = vector("list", R),
    w_store    = vector("list", R),
    b_store    = vector("list", R),
    t_store    = vector("list", R),
    rmse_store = vector("list", R),
    cor_store  = numeric(R),
    Imp_store  = vector("list", R)
  )
  
})
names(results) <- methods

# Store only first 10 datasets globally (not per method)
full_store <- vector("list", 10)
miss_store <- vector("list", 10)

# ============================================================
# TIMING SETUP
# ============================================================

n_cores <- 5   # your machine
cat("Using", n_cores, "cores (manual setting)\n\n")

sim_times    <- numeric(R)
method_times <- matrix(NA, nrow = R, ncol = length(methods))
colnames(method_times) <- methods

global_start <- Sys.time()
# ============================================================
# BACKUP SETUP
# ============================================================

backup_file <- "FirstRealTestn20_backup.rds"

# ============================================================
# SIMULATION LOOP
# ============================================================

for (r in 1:R) {
  
  sim_start <- Sys.time()
  
  cat("====================================================\n")
  cat("Simulation:", r, "of", R, "\n")
  cat("====================================================\n")
  
  # -----------------------------
  # 1) Generate DGP
  # -----------------------------
  
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
  # Derived terms on the FULL data (needed for RMSE + truth comparisons)
  full_df$X3_sq <- full_df$X3^2
  full_df$X1_X2 <- full_df$X1 * full_df$X2
  full_df$X8_X9 <- full_df$X8 * full_df$X9
  # these are just interactions!
  
  
  # -----------------------------
  # 2) MAR Missingness
  # -----------------------------
  
  # miss_df <- full_df
  # p_miss  <- logit_p(0.5 * miss_df$X9 - 0.5 * miss_df$X10)
  # 
  # for (v in c("Y", paste0("X", 1:8))) {
  #   miss_indicator <- rbinom(n, 1, p_miss)
  #   miss_df[miss_indicator == 1, v] <- NA
  # }
  # 
  # if (r <= 10) {
  #   full_store[[r]] <- full_df
  #   miss_store[[r]] <- miss_df
  # }
  # 2) Missingness (MAR vs MCAR)
  miss_df <- full_df
  
  if (missing_mech == "MAR") {
    
    # 50% miss p_miss <- logit_p(0.5 * miss_df$X9 - 0.5 * miss_df$X10)
    # 0.85 notwendig um auf 30% missingh zu kommen
    linpred <- -0.85 + 0.5 * miss_df$X9 - 0.5 * miss_df$X10 
    
    # 10% Missing
    #linpred <- -2.20 + 0.5 * miss_df$X9 - 0.5 * miss_df$X10
    p_miss <- logit_p(linpred)
    #logit^(−1)(−0.85)≈0.30
    
    for (v in c("Y", paste0("X", 1:8))) {
      miss_indicator <- rbinom(n, 1, p_miss)
      miss_df[miss_indicator == 1, v] <- NA
    }
    
  } else if (missing_mech == "MCAR") {
    
    for (v in c("Y", paste0("X", 1:8))) {
      miss_indicator <- rbinom(n, 1, pi_mcar)
      miss_df[miss_indicator == 1, v] <- NA
    }
    
  } else {
    stop("missing_mech must be 'MAR' or 'MCAR'")
  }
  # -------------------------------------------------
  # Recompute derived terms AFTER missingness
  # (so NA structure is correct and no leakage occurs)
  # -------------------------------------------------
  miss_df$X3_sq <- miss_df$X3^2
  miss_df$X1_X2 <- miss_df$X1 * miss_df$X2
  miss_df$X8_X9 <- miss_df$X8 * miss_df$X9
  
  # Store datasets for diagnostics
  if (r <= 10) {
    full_store[[r]] <- full_df
    miss_store[[r]] <- miss_df
  }

  # ============================================================
  # 3) MULTI METHOD IMPUTATION
  # ============================================================
  
  for (method in methods) {
    
    method_start <- Sys.time()
    cat("  Method:", method, "\n")
    
    if (method == "rf_ranger") {
      
      Imp_obj <- miceRanger(
        miss_df,
        m = m_val,
        maxit = 20,
        num.trees = 200,
        num.threads = 5,
        verbose = FALSE
      )
      
      imputed_list <- completeData(Imp_obj)
      
    } else if (method == "pmm") {
      
      Imp_obj <- mice(
        miss_df,
        m = m_val,
        method = "pmm",
        maxit = 20,
        printFlag = FALSE
      )
      
      imputed_list <- complete(Imp_obj, "all")
      
    } else if (method == "cart") {
      
      Imp_obj <- mice(
        miss_df,
        m = m_val,
        method = "cart",
        maxit = 20,
        printFlag = FALSE
      )
      
      imputed_list <- complete(Imp_obj, "all")
    }
    
    results[[method]]$Imp_store[[r]] <- Imp_obj
    
    # -------------------------------------------------
    # Diagnostics
    # -------------------------------------------------
    
    imp1 <- as.data.frame(imputed_list[[1]])
    
    rmse_vals <- sapply(colnames(miss_df), function(v) {
      idx <- which(is.na(miss_df[[v]]))
      if (length(idx) == 0) return(NA_real_)
      sqrt(mean((imp1[[v]][idx] - full_df[[v]][idx])^2))
    })
    
    cor_error <- sqrt(sum((cor(full_df) - cor(imp1))^2))
    
    results[[method]]$rmse_store[[r]] <- rmse_vals
    results[[method]]$cor_store[r]    <- cor_error
    
    # -------------------------------------------------
    # Analysis model + pooling
    # -------------------------------------------------
    
    yfit_list <- lapply(imputed_list, function(dat)
      lm(form_true, data = dat)
    )
    
    class(yfit_list) <- "mira"
    pooled_obj  <- mice::pool(yfit_list)
    pooled_table <- pooled_obj$pooled
    
    pooled_table <- pooled_table[match(names(true_beta),
                                       pooled_table$term), ]
    
    # Store pooled results
    results[[method]]$coef_store[[r]] <- pooled_table$estimate
    results[[method]]$se_store[[r]]   <- sqrt(pooled_table$t)
    results[[method]]$df_store[[r]]   <- pooled_table$df
    results[[method]]$w_store[[r]]    <- pooled_table$ubar
    results[[method]]$b_store[[r]]    <- pooled_table$b
    results[[method]]$t_store[[r]]    <- pooled_table$t
    
    # ---- Method timing ----
    method_end  <- Sys.time()
    method_time <- as.numeric(difftime(method_end,
                                       method_start,
                                       units = "secs"))
    
    method_times[r, method] <- method_time
    cat("    Time (sec):", round(method_time, 2), "\n")
  }
  
  # ============================================================
  # AUTO BACKUP EVERY 10 SIMS
  # ============================================================
  
  if (r %% 10 == 0) {
    
    backup_results <- list(
      results          = results,
      beta             = beta,
      true_beta        = true_beta,
      Sigma            = Sigma,
      form_true        = deparse(form_true),
      n                = n,
      p                = p,
      R_completed      = r,
      m_val            = m_val,
      missing_mechanism= missing_mech,
      seed             = 123,
      full_store       = full_store,
      miss_store       = miss_store,
      sim_times        = sim_times,
      method_times     = method_times,
      date             = Sys.time()
    )
    
    saveRDS(backup_results, backup_file)
    cat(">>> Backup saved at simulation", r, "\n\n")
  }
  
  # ---- Simulation timing ----
  sim_end  <- Sys.time()
  sim_time <- as.numeric(difftime(sim_end,
                                  sim_start,
                                  units = "secs"))
  
  sim_times[r] <- sim_time
  
  elapsed   <- as.numeric(difftime(sim_end,
                                   global_start,
                                   units = "secs"))
  
  avg_time  <- mean(sim_times[1:r])
  est_total <- avg_time * R
  
  cat("Simulation", r, "took", round(sim_time, 2), "sec\n")
  cat("Elapsed:", round(elapsed/60, 2), "min\n")
  cat("Estimated total runtime:",
      round(est_total/60, 2), "min\n\n")
}


cat("\n====================================================\n")
cat("FINAL TIMING SUMMARY\n")
cat("====================================================\n")

cat("Average simulation time (sec):",
    round(mean(sim_times),2), "\n")

cat("Total runtime (min):",
    round(sum(sim_times)/60,2), "\n\n")

cat("Average method times (sec):\n")
print(round(colMeans(method_times),2))


# ============================================================
# EVALUATION PER METHOD
# ============================================================

summary_results <- list()

for (method in methods) {
  
  coef_mat <- do.call(rbind, results[[method]]$coef_store)
  se_mat   <- do.call(rbind, results[[method]]$se_store)
  
  colnames(coef_mat) <- names(true_beta)
  colnames(se_mat)   <- names(true_beta)
  
  bias <- colMeans(coef_mat) - true_beta
  emp_var <- apply(coef_mat, 2, var)
  mean_rubin_var <- colMeans(se_mat^2)
  
  df_mat  <- do.call(rbind, results[[method]]$df_store)
  crit_mat <- qt(0.975, df = df_mat)
  
  tb_mat <- matrix(true_beta, nrow = R, ncol = length(true_beta), byrow = TRUE)
  
  coverage <- colMeans(
    (coef_mat - crit_mat * se_mat <= tb_mat) &
      (coef_mat + crit_mat * se_mat >= tb_mat)
  )
  
  rmse_mat <- do.call(rbind, results[[method]]$rmse_store)
  mean_rmse <- colMeans(rmse_mat, na.rm = TRUE)
  mean_cor_error <- mean(results[[method]]$cor_store)
  
  summary_results[[method]] <- list(
    coef_mat       = coef_mat,
    se_mat         = se_mat,
    bias           = bias,
    emp_var        = emp_var,
    mean_rubin_var = mean_rubin_var,
    coverage       = coverage,
    mean_rmse      = mean_rmse,
    mean_cor_error = mean_cor_error, 
    full_store = full_store,
    miss_store = miss_store
  )
}


# ============================================================
# SAVE
# ============================================================

final_results <- list(
  results         = results,
  summary_results = summary_results,
  beta            = beta,
  true_beta       = true_beta,
  Sigma           = Sigma,
  form_true       = deparse(form_true),
  n               = n,
  p               = p,
  R               = R,
  m_val           = m_val,
  missing_mechanism = "MAR: logit(linpred), linpred = -0.85 + 0.5*X9 - 0.5*X10; derived terms recomputed after missingness",

  seed            = 123,
  full_store      = full_store,   # <-- ADD THIS
  miss_store      = miss_store,   # <-- ADD THIS
  sessionInfo     = sessionInfo(),
  date            = Sys.time()
)

saveRDS(final_results, "FirstRealTestn20.rds")

cat("\nMulti-method simulation complete and saved.\n")
