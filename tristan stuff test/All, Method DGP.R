# ============================================================
# BURGETTE STYLE SIMULATION — MULTI METHOD VERSION
# ============================================================

pacman::p_load(miceRanger, tidyverse, MASS, mice)

set.seed(123)

# -----------------------------
# Parameters
# -----------------------------
n     <- 1000
p     <- 10
R     <- 1
m_val <- 20

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

options(ranger.num.threads = 5)

# ============================================================
# METHODS + STORAGE
# ============================================================

methods <- c("rf_ranger", "pmm", "cart") #, "rf_mice"

results <- lapply(methods, function(x) {
  list(
    coef_store = vector("list", R),
    se_store   = vector("list", R),
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

library(parallel)

n_cores <- 5   # your machine
cat("Using", n_cores, "cores (manual setting)\n\n")

sim_times    <- numeric(R)
method_times <- matrix(NA, nrow = R, ncol = length(methods))
colnames(method_times) <- methods

global_start <- Sys.time()

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
  
  # -----------------------------
  # 2) MAR Missingness
  # -----------------------------
  
  miss_df <- full_df
  p_miss  <- logit_p(0.5 * miss_df$X9 - 0.5 * miss_df$X10)
  
  for (v in c("Y", paste0("X", 1:8))) {
    miss_indicator <- rbinom(n, 1, p_miss)
    miss_df[miss_indicator == 1, v] <- NA
  }
  
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
    }
    
    # else if (method == "rf_mice") {
    #   
    #   Imp_obj <- mice(
    #     miss_df,
    #     m = m_val,
    #     method = "rf",
    #     maxit = 10,
    #     ntree = 200,
    #     printFlag = FALSE
    #   )
    #   
    #   imputed_list <- complete(Imp_obj, "all")
    # }
    
    else if (method == "pmm") {
      
      Imp_obj <- mice(
        miss_df,
        m = m_val,
        method = "pmm",
        maxit = 20,
        printFlag = FALSE
      )
      
      imputed_list <- complete(Imp_obj, "all")
    }
    
    else if (method == "cart") {
      
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
    
    # Diagnostics
    imp1 <- as.data.frame(imputed_list[[1]])
    
    rmse_vals <- sapply(colnames(miss_df), function(v) {
      idx <- which(is.na(miss_df[[v]]))
      if (length(idx) == 0) return(NA_real_)
      sqrt(mean((imp1[[v]][idx] - full_df[[v]][idx])^2))
    })
    
    cor_error <- sqrt(sum((cor(full_df) - cor(imp1))^2))
    
    results[[method]]$rmse_store[[r]] <- rmse_vals
    results[[method]]$cor_store[r]    <- cor_error
    
    # Analysis model
    yfit_list <- lapply(imputed_list, function(dat)
      lm(form_true, data = dat)
    )
    
    class(yfit_list) <- "mira"
    pooled_obj <- mice::pool(yfit_list)
    pooled_table <- pooled_obj$pooled
    
    results[[method]]$coef_store[[r]] <- pooled_table$estimate
    results[[method]]$se_store[[r]]   <- sqrt(pooled_table$t)
    
    # ---- Method timing ----
    method_end <- Sys.time()
    method_time <- as.numeric(difftime(method_end, method_start, units = "secs"))
    
    method_times[r, method] <- method_time
    
    cat("    Time (sec):", round(method_time,2), "\n")
  }
  
  # ---- Simulation timing ----
  sim_end <- Sys.time()
  sim_time <- as.numeric(difftime(sim_end, sim_start, units = "secs"))
  sim_times[r] <- sim_time
  
  elapsed <- as.numeric(difftime(sim_end, global_start, units = "secs"))
  avg_time <- mean(sim_times[1:r])
  est_total <- avg_time * R
  
  cat("Simulation", r, "took", round(sim_time,2), "sec\n")
  cat("Elapsed:", round(elapsed/60,2), "min\n")
  cat("Estimated total runtime:", round(est_total/60,2), "min\n\n")
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
  
  coverage <- colMeans(
    (coef_mat - 1.96 * se_mat <= true_beta) &
      (coef_mat + 1.96 * se_mat >= true_beta)
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
  missing_mechanism = "MAR: logit(0.5*X9 - 0.5*X10)",
  seed            = 123,
  full_store      = full_store,   # <-- ADD THIS
  miss_store      = miss_store,   # <-- ADD THIS
  sessionInfo     = sessionInfo(),
  date            = Sys.time()
)

saveRDS(final_results, "test_burgette_full_sim_multi_method.rds")

cat("\nMulti-method simulation complete and saved.\n")


#############
#Load
burgette_results <- readRDS("test_burgette_full_sim_multi_method.rds")

results         <- burgette_results$results
summary_results <- burgette_results$summary_results
true_beta       <- burgette_results$true_beta
m_val           <- burgette_results$m_val
methods         <- names(summary_results)
full_store      <- burgette_results$full_store
miss_store      <- burgette_results$miss_store
# ============================================================
# RECOVER DATA FOR HTML / EXTRA ANALYSIS
# ============================================================

r <- 1                      # choose simulation
method_html <- methods[1]   # choose method

# Recover full data
full_df <- full_store[[r]]

# Recover imputed list
Imp_obj <- results[[method_html]]$Imp_store[[r]]

if (method_html == "rf_ranger") {
  imputed_list <- completeData(Imp_obj)
} else {
  imputed_list <- complete(Imp_obj, "all")
}

# ============================================================
# PRINT RESULTS (HTML STYLE OUTPUT)
# ============================================================

cat("\n=====================================================\n")
cat("MULTI-METHOD SUMMARY\n")
cat("=====================================================\n\n")

for (method in methods) {
  
  cat("-----------------------------------------------------\n")
  cat("Method:", method, "\n")
  cat("-----------------------------------------------------\n")
  
  summ <- summary_results[[method]]
  
  cat("\nMean Bias:\n")
  print(round(summ$bias, 4))
  
  cat("\nCoverage (should be ~0.95):\n")
  print(round(summ$coverage, 4))
  
  cat("\nEmpirical Variance:\n")
  print(round(summ$emp_var, 4))
  
  cat("\nMean Rubin Variance:\n")
  print(round(summ$mean_rubin_var, 4))
  
  cat("\nMean RMSE (missing entries only):\n")
  print(round(summ$mean_rmse, 4))
  
  cat("\nMean Correlation Matrix Error:\n")
  print(round(summ$mean_cor_error, 4))
  
  cat("\n\n")
}

# ============================================================
# OPTIONAL: BETWEEN-IMPUTATION VARIANCE (B) SUMMARY
# ============================================================

cat("=====================================================\n")
cat("Between-Imputation Variance (Mean across parameters)\n")
cat("=====================================================\n\n")

for (method in methods) {
  
  coef_mat <- summary_results[[method]]$coef_mat
  Q_bar    <- colMeans(coef_mat)
  
  B <- apply(coef_mat, 2, function(q) {
    mean((q - mean(q))^2)
  })
  
  cat("Method:", method, "\n")
  print(round(mean(B), 5))
  cat("\n")
}

# ============================================================
# COMPARISON TABLE
# ============================================================

comparison_table <- data.frame(
  Method      = methods,
  MeanAbsBias = NA_real_,
  Coverage    = NA_real_,
  RMSE_Y      = NA_real_
)

for (i in seq_along(methods)) {
  
  method <- methods[i]
  summ   <- summary_results[[method]]
  
  # Mean absolute bias across parameters
  comparison_table$MeanAbsBias[i] <-
    mean(abs(summ$bias))
  
  # Mean coverage across parameters
  comparison_table$Coverage[i] <-
    mean(summ$coverage)
  
  # RMSE only for Y
  comparison_table$RMSE_Y[i] <-
    summ$mean_rmse["Y"]
}

# round for clean output
comparison_table <- comparison_table %>%
  mutate(
    MeanAbsBias = round(MeanAbsBias, 4),
    Coverage    = round(Coverage, 4),
    RMSE_Y      = round(RMSE_Y, 4)
  )

print(comparison_table)

# -----------------------------
# HTML STYLE QUANTITIES
# -----------------------------

true_mean  <- mean(full_df$Y)
true_quant <- quantile(full_df$Y, 0.9)

theta_MI <- sapply(imputed_list, function(dat) mean(dat$Y))
quant_MI <- sapply(imputed_list, function(dat) quantile(dat$Y, 0.9))

# Rubin components for mean
Q_bar <- mean(theta_MI)
B     <- var(theta_MI)
U_bar <- mean(sapply(imputed_list, function(dat)
  var(dat$Y) / nrow(dat)
))
Tot   <- U_bar + (1 + 1/m_val) * B
se    <- sqrt(Tot)

CIlow  <- Q_bar - 1.96 * se
CIup   <- Q_bar + 1.96 * se

relBias_mean  <- (Q_bar - true_mean) / true_mean
Coverage_mean <- as.numeric(true_mean >= CIlow &
                              true_mean <= CIup)


# Rubin components for quantile
Q_bar_q <- mean(quant_MI)
B_q     <- var(quant_MI)
U_bar_q <- U_bar   # same approx as HTML
Tot_q   <- U_bar_q + (1 + 1/m_val) * B_q
se_q    <- sqrt(Tot_q)

CIlow_q <- Q_bar_q - 1.96 * se_q
CIup_q  <- Q_bar_q + 1.96 * se_q

relBias_quant  <- (Q_bar_q - true_quant) / true_quant
Coverage_quant <- as.numeric(true_quant >= CIlow_q &
                               true_quant <= CIup_q)


# Tables ------------------------------------------------------------------


table_mean <- data.frame(
  Quantity        = "Mean(Y)",
  True_Value      = true_mean,
  Q_bar           = Q_bar,
  CIlow           = CIlow,
  CIup            = CIup,
  Relative_Bias   = relBias_mean,
  Coverage        = Coverage_mean,
  Between_Var_B   = B,
  Within_Var_Ubar = U_bar,
  Total_Var       = Tot,
  SE              = se
)

table_mean[ , -1] <- round(table_mean[ , -1], 4)

print(table_mean)


table_quant <- data.frame(
  Quantity        = "Quantile(0.9)",
  True_Value      = true_quant,
  Q_bar           = Q_bar_q,
  CIlow           = CIlow_q,
  CIup            = CIup_q,
  Relative_Bias   = relBias_quant,
  Coverage        = Coverage_quant,
  Between_Var_B   = B_q,
  Within_Var_Ubar = U_bar_q,
  Total_Var       = Tot_q,
  SE              = se_q
)

table_quant[ , -1] <- round(table_quant[ , -1], 4)

print(table_quant)

table_combined <- rbind(
  data.frame(
    Quantity = "Mean(Y)",
    True     = true_mean,
    Estimate = Q_bar,
    RelBias  = relBias_mean,
    Coverage = Coverage_mean,
    B        = B,
    Ubar     = U_bar,
    TotalVar = Tot
  ),
  data.frame(
    Quantity = "Quantile(0.9)",
    True     = true_quant,
    Estimate = Q_bar_q,
    RelBias  = relBias_quant,
    Coverage = Coverage_quant,
    B        = B_q,
    Ubar     = U_bar_q,
    TotalVar = Tot_q
  )
)

table_combined[ , -1] <- round(table_combined[ , -1], 4)

print(table_combined)


# Plots Bias Rmse etc -----------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

# -------------------------------------------------
# Prepare plotting data (add variance)
# -------------------------------------------------

plot_df_var <- lapply(methods, function(m) {
  
  summ <- summary_results[[m]]
  
  data.frame(
    Method = m,
    Parameter = names(summ$emp_var),
    EmpiricalVar = summ$emp_var,
    RubinVar = summ$mean_rubin_var
  )
  
}) %>% bind_rows()

# Ensure correct parameter order
plot_df_var$Parameter <- factor(
  plot_df_var$Parameter,
  levels = names(true_beta)
)

# Convert to long format for side-by-side bars
plot_df_var_long <- plot_df_var %>%
  pivot_longer(
    cols = c(EmpiricalVar, RubinVar),
    names_to = "VarianceType",
    values_to = "Variance"
  )

# -------------------------------------------------
# Variance Plot
# -------------------------------------------------

p_var <- ggplot(plot_df_var_long,
                aes(x = Parameter,
                    y = Variance,
                    fill = VarianceType)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8)) +
  facet_wrap(~ Method) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Empirical vs Rubin Variance",
       y = "Variance",
       x = "Parameter",
       fill = "")

# -------------------------------------------------
# Show all three stacked
# -------------------------------------------------
p_var
# -------------------------------------------------
# Prepare Bias + Coverage plotting data
# -------------------------------------------------

plot_df <- lapply(methods, function(m) {
  
  summ <- summary_results[[m]]
  
  data.frame(
    Method = m,
    Parameter = names(summ$bias),
    Bias = summ$bias,
    Coverage = summ$coverage
  )
  
}) %>% bind_rows()

plot_df$Parameter <- factor(
  plot_df$Parameter,
  levels = names(true_beta)
)
p_bias <- ggplot(plot_df,
                 aes(x = Parameter,
                     y = Bias,
                     fill = Method)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Bias by Parameter",
       y = "Bias",
       x = "Parameter")

p_cov <- ggplot(plot_df,
                aes(x = Parameter,
                    y = Coverage,
                    fill = Method)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0.95,
             linetype = "dashed",
             color = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Coverage by Parameter",
       y = "Coverage",
       x = "Parameter")

grid.arrange(p_bias, p_cov, ncol = 1)



# Plots first Dataset ------------------------------------------------------

library(gridExtra)
library(grid)

r <- 1
methods_diag <- methods

conv_plots <- list()

for (m in methods_diag) {
  
  Imp_obj <- results[[m]]$Imp_store[[r]]
  
  if (m == "rf_ranger") {
    
    p <- plotVarConvergence(Imp_obj, vars = "Y") +
      ggtitle(paste("Convergence -", m))
    
    conv_plots[[m]] <- p
    
  } else {
    
    p <- grid.grabExpr(
      print(plot(Imp_obj), newpage = FALSE)
    )
    
    p <- arrangeGrob(p,
                     top = textGrob(paste("Convergence -", m),
                                    gp = gpar(fontsize = 14, fontface = "bold")))
    
    conv_plots[[m]] <- p
  }
}

do.call(grid.arrange, c(conv_plots, ncol = 2))

# Print individually
for (m in names(conv_plots)) {
  grid.newpage()
  grid.draw(conv_plots[[m]])
}

library(mice)

plot(Imp_obj, c("Y"))  # default

# Better:
stripplot(Imp_obj, Y ~ .imp)

windows(width = 14, height = 10)

plot(
  Imp_obj,
  c("Y"),
  layout = c(1, 1),
  col = 1:20,
  lwd = 2
)



# ------------------------------------------------------------
# 2) DENSITY PLOTS (2x2)
# ------------------------------------------------------------
dens_plots <- list()

for (m in methods_diag) {
  
  Imp_obj <- results[[m]]$Imp_store[[r]]
  
  if (m == "rf_ranger") {
    
    p <- plotDistributions(Imp_obj, vars = "Y") +
      ggtitle(paste("Density -", m))
    
    dens_plots[[m]] <- p
    
  } else {
    
    p <- grid.grabExpr(
      print(densityplot(Imp_obj, ~ Y), newpage = FALSE)
    )
    
    p <- arrangeGrob(p,
                     top = textGrob(paste("Density -", m),
                                    gp = gpar(fontsize = 14, fontface = "bold")))
    
    dens_plots[[m]] <- p
  }
}

do.call(grid.arrange, c(dens_plots, ncol = 2))



# ------------------------------------------------------------
# 3) STRIPPLOT / IMPUTATION SPREAD (2x2)
# ------------------------------------------------------------
strip_plots <- list()

for (m in methods_diag) {
  
  Imp_obj <- results[[m]]$Imp_store[[r]]
  
  if (m == "rf_ranger") {
    
    p <- plotImputationVariance(Imp_obj, vars = "Y") +
      ggtitle(paste("Imputation Spread -", m))
    
    strip_plots[[m]] <- p
    
  } else {
    
    p <- grid.grabExpr(
      print(stripplot(Imp_obj, Y ~ .imp,
                      pch = 20, cex = 0.6),
            newpage = FALSE)
    )
    
    p <- arrangeGrob(p,
                     top = textGrob(paste("Imputation Spread -", m),
                                    gp = gpar(fontsize = 14, fontface = "bold")))
    
    strip_plots[[m]] <- p
  }
}

do.call(grid.arrange, c(strip_plots, ncol = 2))


# ====================================================
#   Simulation: 1 of 3 
# ====================================================
#   Method: rf_ranger 
# Time (sec): 77.13 
# Method: pmm 
# Time (sec): 5.62 
# Method: cart 
# Time (sec): 40.46 
# Simulation 1 took 123.21 sec
# Elapsed: 2.1 min
# Estimated total runtime: 6.16 min
# 
# ====================================================
#   Simulation: 2 of 3 
# ====================================================
#   Method: rf_ranger 
# Time (sec): 72.04 
# Method: pmm 
# Time (sec): 5.42 
# Method: cart 
# Time (sec): 38.4 
# Simulation 2 took 115.87 sec
# Elapsed: 4.03 min
# Estimated total runtime: 5.98 min
# 
# ====================================================
#   Simulation: 3 of 3 
# ====================================================
#   Method: rf_ranger 
# Time (sec): 69.31 
# Method: pmm 
# Time (sec): 5.47 
# Method: cart 
# Time (sec): 39.81 
# Simulation 3 took 114.59 sec
# Elapsed: 5.94 min
# Estimated total runtime: 5.89 min

