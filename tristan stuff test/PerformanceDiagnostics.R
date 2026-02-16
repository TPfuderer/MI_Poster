
#####################################################
# DGP END 
#########################################################

#############
#Load
burgette_results <- readRDS("FirstRealTestn20.rds")

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

############################################################
# EXTENDED PERFORMANCE DIAGNOSTICS (FIXED)
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

performance_table <- data.frame(
  Method = methods,
  RMSE_Y = NA,
  MeanDiff_Y = NA,
  VarRatio_Y = NA,
  Correlation_Error = NA,
  RelBias_Param = NA,
  SD_vs_SE_Ratio = NA,
  KS_Distance = NA,
  Cov_Frobenius = NA
)

# 🔥 FIX 1: Initialize BEFORE loop
plot_storage_perf <- list()

for (i in seq_along(methods)) {
  
  method <- methods[i]
  summ   <- summary_results[[method]]
  
  Imp_obj <- results[[method]]$Imp_store[[1]]
  
  if (method == "rf_ranger") {
    imputed_list <- completeData(Imp_obj)
  } else {
    imputed_list <- complete(Imp_obj, "all")
  }
  
  full_df <- full_store[[1]]
  full_y  <- full_df$Y
  
  imp_y <- rowMeans(sapply(imputed_list, function(dat) dat$Y))
  
  # --------------------------------------------------
  # Distribution metrics
  # --------------------------------------------------
  
  rmse_y    <- summ$mean_rmse["Y"]
  mean_diff <- mean(imp_y) - mean(full_y)
  var_ratio <- var(imp_y) / var(full_y)
  cor_error <- summ$mean_cor_error
  
  # 🔥 Safer KS (optional but stable)
  ks_dist <- max(abs(ecdf(imp_y)(full_y) - ecdf(full_y)(full_y)))
  
  # 🔥 FIX 2: Proper covariance averaging
  cov_full <- cov(full_df)
  cov_imp_list <- lapply(imputed_list, cov)
  cov_imp  <- Reduce("+", cov_imp_list) / length(cov_imp_list)
  frob_norm <- sqrt(sum((cov_full - cov_imp)^2))
  
  # --------------------------------------------------
  # Parameter-level metrics
  # --------------------------------------------------
  
  coef_mat <- do.call(rbind, results[[method]]$coef_store)
  se_mat   <- do.call(rbind, results[[method]]$se_store)
  
  nonzero_index <- which(true_beta != 0)
  
  rel_bias <- mean(
    abs((colMeans(coef_mat)[nonzero_index] -
           true_beta[nonzero_index]) /
          true_beta[nonzero_index])
  )
  
  emp_sd   <- apply(coef_mat, 2, sd)
  mean_se  <- colMeans(se_mat)
  sd_se_ratio <- mean(mean_se / emp_sd)
  
  # --------------------------------------------------
  # Store table
  # --------------------------------------------------
  
  performance_table$RMSE_Y[i] <- rmse_y
  performance_table$MeanDiff_Y[i] <- mean_diff
  performance_table$VarRatio_Y[i] <- var_ratio
  performance_table$Correlation_Error[i] <- cor_error
  performance_table$RelBias_Param[i] <- rel_bias
  performance_table$SD_vs_SE_Ratio[i] <- sd_se_ratio
  performance_table$KS_Distance[i] <- ks_dist
  performance_table$Cov_Frobenius[i] <- frob_norm
  
  # 🔥 FIX 3: Store for plots
  plot_storage_perf[[method]] <- data.frame(
    Method = method,
    RMSE = rmse_y,
    VarRatio = var_ratio,
    CorError = cor_error
  )
}

performance_table[,-1] <- round(performance_table[,-1], 3)

cat("\n================ EXTENDED PERFORMANCE SUMMARY =================\n")
print(performance_table)
cat("===============================================================\n\n")

# ==========================================================
# PLOTS
# ==========================================================

plot_df_perf <- dplyr::bind_rows(plot_storage_perf)

p_rmse <- ggplot(plot_df_perf,
                 aes(x=Method, y=RMSE, fill=Method)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(title="RMSE (Missing Y only)", y="RMSE")

p_varratio <- ggplot(plot_df_perf,
                     aes(x=Method, y=VarRatio, fill=Method)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=1, linetype="dashed") +
  theme_minimal() +
  labs(title="Variance Ratio (Imputed / True)", y="Ratio")

p_cor <- ggplot(plot_df_perf,
                aes(x=Method, y=CorError, fill=Method)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(title="Correlation Structure Error", y="Error")

gridExtra::grid.arrange(p_rmse, p_varratio, p_cor, ncol=1)

############################################################
# RUBIN DIAGNOSTICS — SINGLE DATASET (ALL METHODS)
############################################################

r <- 1   # chosen dataset

rubin_compare <- list()

for (method in methods) {
  
  Imp_obj <- results[[method]]$Imp_store[[r]]
  
  if (method == "rf_ranger") {
    imputed_list <- completeData(Imp_obj)
  } else {
    imputed_list <- complete(Imp_obj, "all")
  }
  
  m <- length(imputed_list)
  
  B_vec  <- results[[method]]$b_store[[r]]
  T_vec  <- results[[method]]$t_store[[r]]
  SE_vec <- results[[method]]$se_store[[r]]
  
  Ubar_vec <- SE_vec^2
  
  RIV_vec <- ((1 + 1/m) * B_vec) / Ubar_vec
  FMI_vec <- ((1 + 1/m) * B_vec) / T_vec
  df_vec  <- (m - 1) * (1 + 1/RIV_vec)^2
  
  MC_error_vec  <- sqrt(B_vec / m)
  MC_ratio_vec  <- MC_error_vec / sqrt(T_vec)
  
  B_over_U_vec  <- B_vec / Ubar_vec
  
  rubin_compare[[method]] <- data.frame(
    Method = method,
    Parameter = names(true_beta),
    B_over_U = B_over_U_vec,
    FMI = FMI_vec,
    RIV = RIV_vec,
    df = df_vec,
    MC_Error_Ratio = MC_ratio_vec
  )
}

rubin_compare_df <- dplyr::bind_rows(rubin_compare)

rubin_compare_df <- rubin_compare_df |>
  dplyr::mutate(
    B_over_U = round(B_over_U,3),
    FMI = round(FMI,3),
    RIV = round(RIV,3),
    df = round(df,1),
    MC_Error_Ratio = round(MC_Error_Ratio,3)
  )

cat("\n================ RUBIN SINGLE-DATASET COMPARISON =================\n")
print(rubin_compare_df)
cat("===================================================================\n\n")

############################################################
# ML-STYLE METRICS (Not MI-appropriate)
############################################################

ml_metrics <- data.frame(
  Method = methods,
  R_squared = NA,
  MSE_missing = NA
)

for (i in seq_along(methods)) {
  
  method <- methods[i]
  
  Imp_obj <- results[[method]]$Imp_store[[1]]
  
  if (method == "rf_ranger") {
    imputed_list <- completeData(Imp_obj)
  } else {
    imputed_list <- complete(Imp_obj, "all")
  }
  
  full_df <- full_store[[1]]
  miss_df <- miss_store[[1]]
  
  missing_index <- which(is.na(miss_df$Y))
  
  imp_y_missing <- rowMeans(
    sapply(imputed_list, function(dat) dat$Y[missing_index])
  )
  
  true_missing <- full_df$Y[missing_index]
  
  mse <- mean((imp_y_missing - true_missing)^2)
  
  r2  <- cor(imp_y_missing, true_missing)^2
  
  ml_metrics$R_squared[i] <- r2
  ml_metrics$MSE_missing[i] <- mse
}

ml_metrics[,-1] <- round(ml_metrics[,-1],3)

cat("\n================ ML-STYLE METRICS (FOR FUN) =================\n")
print(ml_metrics)
cat("==============================================================\n\n")

