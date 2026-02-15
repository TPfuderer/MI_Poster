# ============================================================
# BURGETTE STYLE SIMULATION — MULTI METHOD VERSION (PARALLEL)
# ============================================================

pacman::p_load(miceRanger, tidyverse, MASS, mice, parallel)

set.seed(123)

# -----------------------------
# Parameters
# -----------------------------
n     <- 1000
p     <- 10
R     <- 3
m_val <- 20

beta <- c(
  b0 = 0, b1 = 0.5, b2 = 0.5, b3 = 0.5,
  b4 = 0.5, b5 = 0.5, b6 = 1,
  b7 = 1, b8 = 1
)

Sigma <- diag(p)
Sigma[1:4, 1:4]   <- 0.5; diag(Sigma[1:4, 1:4])   <- 1
Sigma[5:10, 5:10] <- 0.3; diag(Sigma[5:10, 5:10]) <- 1

logit_p <- function(z) 1 / (1 + exp(-z))

form_true <- Y ~ X1 + X2 + X3 + X8 + X9 + I(X3^2) + X1:X2 + X8:X9

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

methods <- c("rf_ranger", "pmm", "cart")

options(ranger.num.threads = 1)  # IMPORTANT inside parallel!

# ============================================================
# PARALLEL SETUP
# ============================================================

n_cores <- 5
cl <- makeCluster(n_cores)
clusterSetRNGStream(cl, 123)

clusterExport(cl, ls())

cat("Running", R, "simulations on", n_cores, "cores...\n\n")

global_start <- Sys.time()

# ============================================================
# PARALLEL SIMULATION
# ============================================================

sim_results <- parLapply(cl, 1:R, function(r) {
  
  cat("Simulation", r, "started\n")
  
  library(miceRanger)
  library(MASS)
  library(mice)
  
  sim_start <- Sys.time()
  
  # -----------------------------
  # 1) Generate Data
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
  
  # -----------------------------
  # 2) Missingness
  # -----------------------------
  miss_df <- full_df
  p_miss <- logit_p(0.5 * miss_df$X9 - 0.5 * miss_df$X10)
  
  for (v in c("Y", paste0("X", 1:8))) {
    miss_indicator <- rbinom(n, 1, p_miss)
    miss_df[miss_indicator == 1, v] <- NA
  }
  
  method_store <- list()
  
  # -----------------------------
  # 3) Methods
  # -----------------------------
  for (method in methods) {
    
    method_start <- Sys.time()
    
    if (method == "rf_ranger") {
      Imp_obj <- miceRanger(
        miss_df,
        m = m_val,
        maxit = 10,
        num.trees = 200,
        num.threads = 1,
        verbose = FALSE
      )
      imputed_list <- completeData(Imp_obj)
    }
    
    else if (method == "pmm") {
      Imp_obj <- mice(
        miss_df,
        m = m_val,
        method = "pmm",
        maxit = 10,
        printFlag = FALSE
      )
      imputed_list <- complete(Imp_obj, "all")
    }
    
    else if (method == "cart") {
      Imp_obj <- mice(
        miss_df,
        m = m_val,
        method = "cart",
        maxit = 10,
        printFlag = FALSE
      )
      imputed_list <- complete(Imp_obj, "all")
    }
    
    # Diagnostics
    imp1 <- as.data.frame(imputed_list[[1]])
    
    rmse_vals <- sapply(colnames(miss_df), function(v) {
      idx <- which(is.na(miss_df[[v]]))
      if (length(idx) == 0) return(NA_real_)
      sqrt(mean((imp1[[v]][idx] - full_df[[v]][idx])^2))
    })
    
    cor_error <- sqrt(sum((cor(full_df) - cor(imp1))^2))
    
    yfit_list <- lapply(imputed_list, function(dat)
      lm(form_true, data = dat)
    )
    
    class(yfit_list) <- "mira"
    pooled_obj <- mice::pool(yfit_list)
    pooled_table <- pooled_obj$pooled
    
    method_time <- as.numeric(difftime(Sys.time(), method_start, units = "secs"))
    
    method_store[[method]] <- list(
      coef = pooled_table$estimate,
      se   = sqrt(pooled_table$t),
      rmse = rmse_vals,
      cor  = cor_error,
      time = method_time
    )
  }
  
  sim_time <- as.numeric(difftime(Sys.time(), sim_start, units = "secs"))
  
  cat("Simulation", r, "finished in", round(sim_time,2), "sec\n")
  
  return(list(
    methods = method_store,
    sim_time = sim_time
  ))
})

stopCluster(cl)

cat("\nAll simulations completed.\n")
cat("Total runtime (min):",
    round(as.numeric(difftime(Sys.time(), global_start, units="secs"))/60,2),
    "\n")

# ============================================================
# SAVE
# ============================================================

final_results <- list(
  sim_results = sim_results,
  beta        = beta,
  true_beta   = true_beta,
  Sigma       = Sigma,
  form_true   = deparse(form_true),
  n           = n,
  p           = p,
  R           = R,
  m_val       = m_val,
  seed        = 123,
  sessionInfo = sessionInfo(),
  date        = Sys.time()
)

saveRDS(final_results, "burgette_parallel_sim.rds")

cat("Results saved.\n")
