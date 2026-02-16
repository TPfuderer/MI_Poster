
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
# EXTRA MI DIAGNOSTICS (MCAR vs MAR check)
# ============================================================

for (method in methods) {
  
  coef_mat <- do.call(rbind, results[[method]]$coef_store)
  se_mat   <- do.call(rbind, results[[method]]$se_store)
  df_mat   <- do.call(rbind, results[[method]]$df_store)
  b_mat    <- do.call(rbind, results[[method]]$b_store)
  t_mat    <- do.call(rbind, results[[method]]$t_store)
  
  colnames(coef_mat) <- names(true_beta)
  colnames(se_mat)   <- names(true_beta)
  
  # t critical values
  crit_mat <- qt(0.975, df = df_mat)
  
  tb_mat <- matrix(true_beta, nrow = R,
                   ncol = length(true_beta),
                   byrow = TRUE)
  
  coverage_t <- colMeans(
    (coef_mat - crit_mat * se_mat <= tb_mat) &
      (coef_mat + crit_mat * se_mat >= tb_mat)
  )
  
  # Fraction of Missing Information
  FMI_mat <- ((1 + 1/m_val) * b_mat) / t_mat
  mean_FMI <- colMeans(FMI_mat, na.rm = TRUE)
  
  # Variance comparison
  emp_var <- apply(coef_mat, 2, var)
  mean_T  <- colMeans(t_mat, na.rm = TRUE)
  
  cat("\n==================== DIAGNOSTICS:", method, "====================\n")
  cat("Missingness mech:", missing_mech, "\n")
  cat("Missing rate Y:",
      round(mean(is.na(miss_store[[1]]$Y)), 3), "\n")
  cat("Mean df:",
      round(mean(df_mat, na.rm=TRUE),2),
      "| Mean t-crit:",
      round(mean(crit_mat, na.rm=TRUE),3), "\n")
  cat("Mean coverage (t):",
      round(mean(coverage_t),3), "\n")
  cat("Mean FMI:",
      round(mean(mean_FMI),3), "\n")
  cat("Mean(mean_T / emp_var):",
      round(mean(mean_T / emp_var, na.rm=TRUE),3), "\n")
  cat("===============================================================\n")
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
    TotalVariance = summ$mean_rubin_var
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
    cols = c(EmpiricalVar, TotalVariance),
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



############################################################
# SIMULATION-LEVEL PROPERNESS DIAGNOSTICS
############################################################

for (method in methods) {
  
  coef_mat <- do.call(rbind, results[[method]]$coef_store)
  se_mat   <- do.call(rbind, results[[method]]$se_store)
  df_mat   <- do.call(rbind, results[[method]]$df_store)
  b_mat    <- do.call(rbind, results[[method]]$b_store)
  t_mat    <- do.call(rbind, results[[method]]$t_store)
  
  R_sim <- nrow(coef_mat)
  
  colnames(coef_mat) <- names(true_beta)
  colnames(se_mat)   <- names(true_beta)
  
  crit_mat <- qt(0.975, df = df_mat)
  
  tb_mat <- matrix(true_beta,
                   nrow = R_sim,
                   ncol = length(true_beta),
                   byrow = TRUE)
  
  coverage_t <- colMeans(
    (coef_mat - crit_mat * se_mat <= tb_mat) &
      (coef_mat + crit_mat * se_mat >= tb_mat)
  )
  
  FMI_mat <- ((1 + 1/m_val) * b_mat) / t_mat
  mean_FMI <- colMeans(FMI_mat, na.rm = TRUE)
  
  emp_var <- apply(coef_mat, 2, var)
  mean_T  <- colMeans(t_mat, na.rm = TRUE)
  
  cat("\n==================== SIM DIAGNOSTICS:", method, "====================\n")
  cat("Mean df:", round(mean(df_mat, na.rm=TRUE),2), "\n")
  cat("Mean coverage:", round(mean(coverage_t),3), "\n")
  cat("Mean FMI:", round(mean(mean_FMI),3), "\n")
  cat("Mean(T / empirical var):",
      round(mean(mean_T / emp_var, na.rm=TRUE),3), "\n")
  cat("===============================================================\n")
}

############################################################
# MODEL-SPECIFIC DIAGNOSTICS (Single r)
############################################################

r <- 1
method <- methods[1]

Imp_obj <- results[[method]]$Imp_store[[r]]

cat("\n==================== MODEL CHECK ====================\n")
cat("Method:", method, "\n")

if (method == "rf_ranger") {
  
  print(plotVarConvergence(Imp_obj, vars = "Y"))
  print(plotDistributions(Imp_obj, vars = "Y"))
  print(plotImputationVariance(Imp_obj, vars = "Y"))
  
} else {
  
  plot(Imp_obj)
  densityplot(Imp_obj, ~ Y)
  stripplot(Imp_obj, Y ~ .imp)
}

cat("=====================================================\n")

############################################################
# RUBIN PROPERNESS DIAGNOSTICS
############################################################

rubin_table <- data.frame(
  Method = methods,
  Mean_Coverage = NA,
  Mean_FMI      = NA,
  Mean_T_over_EmpVar = NA,
  Mean_df       = NA
)

for (i in seq_along(methods)) {
  
  method <- methods[i]
  
  coef_mat <- do.call(rbind, results[[method]]$coef_store)
  df_mat   <- do.call(rbind, results[[method]]$df_store)
  b_mat    <- do.call(rbind, results[[method]]$b_store)
  t_mat    <- do.call(rbind, results[[method]]$t_store)
  
  R_sim <- nrow(coef_mat)
  
  crit_mat <- qt(0.975, df = df_mat)
  
  tb_mat <- matrix(true_beta,
                   nrow = R_sim,
                   ncol = length(true_beta),
                   byrow = TRUE)
  
  se_mat <- sqrt(t_mat)
  
  coverage <- colMeans(
    (coef_mat - crit_mat * se_mat <= tb_mat) &
      (coef_mat + crit_mat * se_mat >= tb_mat)
  )
  
  FMI_mat <- ((1 + 1/m_val) * b_mat) / t_mat
  
  emp_var <- apply(coef_mat, 2, var)
  mean_T  <- colMeans(t_mat)
  
  rubin_table$Mean_Coverage[i] <-
    mean(coverage)
  
  rubin_table$Mean_FMI[i] <-
    mean(FMI_mat, na.rm=TRUE)
  
  rubin_table$Mean_T_over_EmpVar[i] <-
    mean(mean_T / emp_var, na.rm=TRUE)
  
  rubin_table$Mean_df[i] <-
    mean(df_mat, na.rm=TRUE)
}

rubin_table[ , -1] <- round(rubin_table[ , -1], 3)
print(rubin_table)


############################################################
# PERFORMANCE DIAGNOSTICS (Not Rubin)
############################################################

perf_table <- data.frame(
  Method      = methods,
  MeanAbsBias = NA,
  RMSE_Y      = NA,
  CorError    = NA
)

for (i in seq_along(methods)) {
  
  method <- methods[i]
  summ   <- summary_results[[method]]
  
  perf_table$MeanAbsBias[i] <-
    mean(abs(summ$bias))
  
  perf_table$RMSE_Y[i] <-
    summ$mean_rmse["Y"]
  
  perf_table$CorError[i] <-
    summ$mean_cor_error
}

perf_table[ , -1] <- round(perf_table[ , -1], 3)
print(perf_table)

mean_T <- colMeans(t_mat)
mean_T
