# ============================================================
# BURGETTE STYLE SIMULATION — MULTI METHOD VERSION (FIXED)
# ============================================================

pacman::p_load(miceRanger, tidyverse, MASS, mice, foreach, doParallel)

set.seed(123)

missing_mech <- "MAR"   # "MAR" or "MCAR"
# pi_mcar <- 0.20       # target missing prob for MCAR
# just used for test mcar not implemented
# -----------------------------
# Parameters
# -----------------------------
n          <- 1000 # Anzahl der Rows
p          <- 10   # Anzahl der X_i Werte
R          <- 300    # Anzahl der Simulationen
m_val      <- 10   # Anzahl der DFs pro Simulation
iterations <- 25   # The number of iterations to run for each dataset.
n_cores    <- 5    # Mögen die Maschinengeister uns beistehen

# -----------------------------
# True coefficients (clean: X1..X10 + 3 derived terms)
# -----------------------------
beta <- c(
  b0    = 0,
  b1    = 0.5,  # X1
  b2    = 0.5,  # X2
  b3    = 0.5,  # X3
  b4    = 0.5,  # X4
  b5    = 0.5,  # X5
  b6    = 0.5,  # X6
  b7    = 0.5,  # X7
  b8    = 0.5,  # X8
  b9    = 0.5,  # X9
  b10   = 0.5,  # X10
  # Derived terms (Indizes nach Variablenname, nicht nach Nummer)
  b_X3sq  = 0.5,  # X3_sq
  b_X1X2  = 1,    # X1_X2
  b_X4X5  = 1     # X4_X5
)

# -----------------------------
# Covariance design (clean)
# X1..X5 correlated 0.5
# X6..X10 correlated 0.3
# -----------------------------
Sigma <- diag(p)
Sigma[1:5, 1:5]   <- 0.5; diag(Sigma[1:5, 1:5])   <- 1
Sigma[6:10, 6:10] <- 0.3; diag(Sigma[6:10, 6:10]) <- 1

logit_p <- function(z) 1 / (1 + exp(-z))

# -----------------------------
# Substantive (analysis) model
# NOTE: derived terms are computed AFTER imputation (fix)
# Erweiterung B: X6-X10 ins Analysemodell aufgenommen
# Benötigen diese für das Analysis Model
# -----------------------------
form_true <- Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X3_sq + X1_X2 + X4_X5

true_beta <- c(
  "(Intercept)" = unname(beta["b0"]),
  "X1"          = unname(beta["b1"]),
  "X2"          = unname(beta["b2"]),
  "X3"          = unname(beta["b3"]),
  "X4"          = unname(beta["b4"]),
  "X5"          = unname(beta["b5"]),
  "X6"          = unname(beta["b6"]),
  "X7"          = unname(beta["b7"]),
  "X8"          = unname(beta["b8"]),
  "X9"          = unname(beta["b9"]),
  "X10"         = unname(beta["b10"]),
  "X3_sq"       = unname(beta["b_X3sq"]),
  "X1_X2"       = unname(beta["b_X1X2"]),
  "X4_X5"       = unname(beta["b_X4X5"])
)

# ============================================================
# METHODS + STORAGE
# ============================================================

methods <- c("rf_ranger", "pmm", "cart")

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

# Store only first 4 datasets globally
full_store <- vector("list", 4)
miss_store <- vector("list", 4)

# ============================================================
# TIMING SETUP
# ============================================================

sim_times    <- numeric(R)
method_times <- matrix(NA, nrow = R, ncol = length(methods))
colnames(method_times) <- methods

global_start <- Sys.time()

# ============================================================
# Helper: add deterministic derived terms (critical fix)
# ============================================================
add_derived <- function(dat){
  dat <- as.data.frame(dat)
  dat$X3_sq <- dat$X3^2
  dat$X1_X2 <- dat$X1 * dat$X2
  dat$X4_X5 <- dat$X4 * dat$X5
  dat
}

# ============================================================
# PARALLEL CLUSTER SETUP
# num.threads = 1 pro Worker: Parallelität läuft auf Prozessebene
# (n_cores Prozesse parallel, jeder single-threaded)
# nested Parallelitaet wuerde zu Ressourcenkonflikten fuehren
# ============================================================

cl <- makeCluster(n_cores)
registerDoParallel(cl)

cat("Using", n_cores, "cores\n\n")

# ============================================================
# SIMULATION LOOP (PARALLEL)
# ============================================================

sim_results <- foreach(
  r         = 1:R,
  .packages = c("miceRanger", "mice", "MASS", "dplyr")
) %dopar% {

  set.seed(r * 100 + 123)
  iter_start <- Sys.time()

  # -----------------------------
  # 1) Generate DGP
  # -----------------------------
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(X) <- paste0("X", 1:10)

  epsilon <- rnorm(n)

  Y <- beta["b0"] +
    beta["b1"]     * X[, "X1"] +
    beta["b2"]     * X[, "X2"] +
    beta["b3"]     * X[, "X3"] +
    beta["b4"]     * X[, "X4"] +
    beta["b5"]     * X[, "X5"] +
    beta["b6"]     * X[, "X6"] +
    beta["b7"]     * X[, "X7"] +
    beta["b8"]     * X[, "X8"] +
    beta["b9"]     * X[, "X9"] +
    beta["b10"]    * X[, "X10"] +
    beta["b_X3sq"] * X[, "X3"]^2 +
    beta["b_X1X2"] * X[, "X1"] * X[, "X2"] +
    beta["b_X4X5"] * X[, "X4"] * X[, "X5"] +
    epsilon

  full_df <- data.frame(Y, X)

  # squared und Interaktionsterme zum df hinzufügen (derived terms)
  full_df <- add_derived(full_df)

  # -----------------------------
  # 2) Missingness (MAR vs MCAR)
  # IMPORTANT: impute only BASE vars; derived are not imputed (fix)
  # -----------------------------
  miss_df <- full_df

  if (missing_mech == "MAR") {

    # MAR: Fehlwahrscheinlichkeit haengt nur von X9 und X10 ab
    # X9 und X10 koennen selbst keine NAs haben
    linpred <- -0.85 + 0.5 * miss_df$X9 - 0.5 * miss_df$X10
    p_miss  <- logit_p(linpred)

    for (v in c("Y", paste0("X", 1:8))) {
      # jede Reihe erhaelt eine eigene Wahrscheinlichkeit, das ein Eintrag darin
      # zu einem NA wird (weil abhaengig nur von X9 und X10 von der Zeile)
      # In diesem Abschnitt werden fuer jede Spalte die zuvor
      # berechnete Wahrscheinlichkeit (innerhalb einer Zeile)
      # in einer Bernoulliverteilung verwendet fuer die Generierung der NA
      # Fuer die Interaktionsterme werden die nicht berechnet, wuerde auch keinen Sinn machen,
      # wenn in einer Zeile X3 existiert aber X3^2 ist ein NA
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

  # NAs in derived terms ergeben sich automatisch aus NAs in Basisvariablen
  # Diese Spalten werden NICHT imputiert
  miss_df <- add_derived(miss_df)

  # Base-only data used for imputation (critical fix)
  miss_imp  <- miss_df %>% dplyr::select(Y, paste0("X", 1:10))
  full_base <- full_df %>% dplyr::select(Y, paste0("X", 1:10))

  # ============================================================
  # 3) MULTI METHOD IMPUTATION
  # ============================================================

  method_results <- list()
  method_times_r <- numeric(length(methods))
  names(method_times_r) <- methods

  for (method in methods) {

    method_start <- Sys.time()

    if (method == "rf_ranger") {

      Imp_obj <- miceRanger::miceRanger(
        miss_imp,
        m           = m_val,
        maxit       = iterations, # bei jedem gleich
        num.trees   = 200, # Standardmäßig mehr => wegen Performance reduziert
        num.threads = 1,  # 1 Thread pro Worker (Parallelität auf Prozessebene)
        verbose     = FALSE
      )
      imputed_list <- miceRanger::completeData(Imp_obj)

    } else if (method == "pmm") {

      Imp_obj <- mice::mice(
        miss_imp,
        m         = m_val,
        method    = "pmm",
        maxit     = iterations,
        printFlag = FALSE
      )
      imputed_list <- mice::complete(Imp_obj, "all")

    } else if (method == "cart") {

      Imp_obj <- mice::mice(
        miss_imp,
        m         = m_val,
        method    = "cart",
        maxit     = iterations,
        printFlag = FALSE
      )
      imputed_list <- mice::complete(Imp_obj, "all")
    }

    # Add derived terms AFTER completion (identities now always hold)
    imputed_list <- lapply(imputed_list, add_derived)

    # -------------------------------------------------
    # Diagnostics
    # -------------------------------------------------
    imp1 <- as.data.frame(imputed_list[[1]])

    # RMSE computed only for base vars (since only those are truly imputed)
    rmse_vals <- sapply(colnames(miss_imp), function(v) {
      idx <- which(is.na(miss_imp[[v]]))
      if (length(idx) == 0) return(NA_real_)
      sqrt(mean((imp1[[v]][idx] - full_base[[v]][idx])^2))
    })

    # Correlation error on full set including derived terms
    full_for_cor <- add_derived(full_base)
    cor_error <- sqrt(sum((cor(full_for_cor) - cor(imp1))^2))

    # -------------------------------------------------
    # Analysis model + pooling
    # -------------------------------------------------
    yfit_list <- lapply(imputed_list, function(dat) lm(form_true, data = dat))
    class(yfit_list) <- "mira"

    pooled_obj   <- mice::pool(yfit_list)
    pooled_table <- pooled_obj$pooled
    pooled_table <- pooled_table[match(names(true_beta), pooled_table$term), ]

    method_results[[method]] <- list(
      coef = pooled_table$estimate,
      se   = sqrt(pooled_table$t),
      df   = pooled_table$df,
      w    = pooled_table$ubar,
      b    = pooled_table$b,
      t    = pooled_table$t,
      rmse = rmse_vals,
      cor  = cor_error,
      Imp  = Imp_obj
    )

    method_times_r[method] <- as.numeric(difftime(Sys.time(), method_start, units = "secs"))
  }

  iter_time <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))

  # Rueckgabewert jeder parallelen Iteration
  list(
    r              = r,
    method_results = method_results,
    method_times   = method_times_r,
    sim_time       = iter_time,
    full_df        = if (r <= 4) full_df else NULL,
    miss_df        = if (r <= 4) miss_df else NULL
  )
}

stopCluster(cl)

# ============================================================
# ERGEBNISSE ZUSAMMENFUEHREN
# foreach gibt eine Liste zurueck — hier werden die Einzelergebnisse
# in die bestehende results-Struktur uebertragen
# ============================================================

for (res in sim_results) {
  r <- res$r
  for (method in methods) {
    mr <- res$method_results[[method]]
    results[[method]]$coef_store[[r]] <- mr$coef
    results[[method]]$se_store[[r]]   <- mr$se
    results[[method]]$df_store[[r]]   <- mr$df
    results[[method]]$w_store[[r]]    <- mr$w
    results[[method]]$b_store[[r]]    <- mr$b
    results[[method]]$t_store[[r]]    <- mr$t
    results[[method]]$rmse_store[[r]] <- mr$rmse
    results[[method]]$cor_store[r]    <- mr$cor
    results[[method]]$Imp_store[[r]]  <- mr$Imp
  }
  sim_times[res$r]      <- res$sim_time
  method_times[res$r, ] <- res$method_times
  if (!is.null(res$full_df)) full_store[[res$r]] <- res$full_df
  if (!is.null(res$miss_df)) miss_store[[res$r]] <- res$miss_df
}

# ============================================================
# TIMING SUMMARY
# ============================================================

cat("\n====================================================\n")
cat("FINAL TIMING SUMMARY\n")
cat("====================================================\n")



total_elapsed <- as.numeric(difftime(Sys.time(), global_start, units = "mins"))
cat("Using", n_cores, "cores\n")

cat("Total CPU time (min):",     round(sum(sim_times) / 60, 2), "\n")
cat("Wall-clock runtime (min):", round(total_elapsed, 2), "\n")                                   
cat("Total CPU time (min):",     round(sum(sim_times) / 60, 2), "\n\n")    
cat("Average method times (sec):\n")
print(round(colMeans(method_times), 2))

# ============================================================
# EVALUATION PER METHOD
# ============================================================
# Damit wir nachher die Plots schneller ersteellen können.
summary_results <- list()

for (method in methods) {

  coef_mat <- do.call(rbind, results[[method]]$coef_store)
  se_mat   <- do.call(rbind, results[[method]]$se_store)

  colnames(coef_mat) <- names(true_beta)
  colnames(se_mat)   <- names(true_beta)

  bias           <- colMeans(coef_mat) - true_beta
  emp_var        <- apply(coef_mat, 2, var)
  mean_rubin_var <- colMeans(se_mat^2)

  df_mat   <- do.call(rbind, results[[method]]$df_store)
  crit_mat <- qt(0.975, df = df_mat)

  tb_mat <- matrix(true_beta, nrow = R, ncol = length(true_beta), byrow = TRUE)

  coverage <- colMeans(
    (coef_mat - crit_mat * se_mat <= tb_mat) &
      (coef_mat + crit_mat * se_mat >= tb_mat)
  )

  rmse_mat       <- do.call(rbind, results[[method]]$rmse_store)
  mean_rmse      <- colMeans(rmse_mat, na.rm = TRUE)
  mean_cor_error <- mean(results[[method]]$cor_store)

  summary_results[[method]] <- list(
    coef_mat       = coef_mat,
    se_mat         = se_mat,
    bias           = bias,
    emp_var        = emp_var,
    mean_rubin_var = mean_rubin_var,
    coverage       = coverage,
    mean_rmse      = mean_rmse,
    mean_cor_error = mean_cor_error
  )
}

# ============================================================
# SAVE
# ============================================================

final_results <- list(
  results           = results,
  summary_results   = summary_results,
  beta              = beta,
  true_beta         = true_beta,
  Sigma             = Sigma,
  form_true         = deparse(form_true),
  n                 = n,
  p                 = p,
  R                 = R,
  m_val             = m_val,
  missing_mechanism = "MAR: logit(linpred), linpred = -0.85 + 0.5*X9 - 0.5*X10; base-only imputation; derived terms computed post-imputation; X6-X10 sind echte Praediktoren (Erweiterung B)",
  seed              = 123,
  full_store        = full_store,
  miss_store        = miss_store,
  sessionInfo       = sessionInfo(),
  date              = Sys.time()
)

saveRDS(final_results, "25_iterations_300sims_uff.rds")

cat("\nMulti-method simulation complete and saved.\n")
