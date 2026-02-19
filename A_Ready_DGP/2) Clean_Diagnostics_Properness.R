
#####################################################
# DGP END 
#########################################################
library(miceRanger)
#############
#Load
burgette_results <- readRDS("Test4.rds")

results         <- burgette_results$results
summary_results <- burgette_results$summary_results
true_beta       <- burgette_results$true_beta
m_val           <- burgette_results$m_val
methods         <- names(summary_results)
full_store      <- burgette_results$full_store
miss_store      <- burgette_results$miss_store

# ============================================================
# RECOVER DATA FOR SINGLE simrun analysis
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
# PROPERNESS DIAGNOSTICS — CLEAN VERSION (ALL SIM RUNS)
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

properness_table <- data.frame(
  Method = methods,
  Mean_Bias = NA,
  Mean_Coverage = NA,
  Mean_T_over_EmpVar = NA
)

plot_storage <- list()

for (i in seq_along(methods)) {
  
  method <- methods[i]
  
  coef_mat <- do.call(rbind, results[[method]]$coef_store)
  t_mat    <- do.call(rbind, results[[method]]$t_store)
  df_mat   <- do.call(rbind, results[[method]]$df_store)
  
  R_sim <- nrow(coef_mat)
  colnames(coef_mat) <- names(true_beta)
  
  # ======================================================
  # 1️⃣ BIAS
  # ======================================================
  
  bias_vec  <- colMeans(coef_mat) - true_beta
  mean_bias <- mean(abs(bias_vec))
  # ======================================================
  # 2️⃣ COVERAGE
  # ======================================================
  
  se_mat   <- sqrt(t_mat)
  crit_mat <- qt(0.975, df = df_mat)
  
  tb_mat <- matrix(true_beta,
                   nrow = R_sim,
                   ncol = length(true_beta),
                   byrow = TRUE)
  
  coverage <- colMeans(
    (coef_mat - crit_mat * se_mat <= tb_mat) &
      (coef_mat + crit_mat * se_mat >= tb_mat)
  )
  
  mean_coverage <- mean(coverage)
  
  # ======================================================
  # 3️⃣ VARIANCE CALIBRATION
  # ======================================================
  
  emp_var <- apply(coef_mat, 2, var)
  mean_T  <- colMeans(t_mat)
  T_ratio <- mean(mean_T / emp_var)
  
  # Store summary
  properness_table$Mean_Bias[i] <- mean_bias
  properness_table$Mean_Coverage[i] <- mean_coverage
  properness_table$Mean_T_over_EmpVar[i] <- T_ratio
  
  # Store parameter-level diagnostics
  plot_storage[[method]] <- data.frame(
    Method = method,
    Parameter = names(true_beta),
    Bias = bias_vec,
    Coverage = coverage,
    EmpVar = emp_var,
    RubinVar = mean_T
  )
}

properness_table[,-1] <- round(properness_table[,-1], 3)

#print(properness_table)

plot_df <- bind_rows(plot_storage)

ratio_df <- plot_df %>%
  mutate(Ratio = RubinVar / EmpVar,
         Calibration = case_when(
           Ratio < 0.95 ~ "Underestimation",
           Ratio > 1.05 ~ "Overestimation",
           TRUE ~ "Well calibrated"
         ))
###################
#Plots
###################


#1️⃣ COVERAGE CHECK
cat("\nCOVERAGE CHECK\n")
cat("Good: ~0.95 | <0.90 = undercoverage | >0.97 overly conservative\n\n")

p_cov <- ggplot(plot_df,
                aes(x = Parameter,
                    y = Coverage,
                    fill = Method)) +
  geom_bar(stat="identity",
           position=position_dodge()) +
  geom_hline(yintercept=0.95,
             linetype="dashed",
             color="red") +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, size = 26, face = "bold"),
    axis.text.y  = element_text(size = 26),
    axis.title   = element_text(size = 30),
    plot.title   = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 28),
    legend.text  = element_text(size = 26)
  ) +
  labs(title = "Coverage by Parameter",
       y = "Coverage",
       x = "Parameter")

print(p_cov)
ggsave("coverage_plot.png",
       plot = p_cov,
       width = 12,
       height = 8,
       dpi = 600)


cat("\n================ COVERAGE =================\n")
cat("Good ≈ 0.95 | <0.90 undercoverage | >0.97 conservative\n\n")

coverage_table <- plot_df %>%
  dplyr::select(Method, Parameter, Coverage) %>%
  pivot_wider(names_from = Method, values_from = Coverage) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

print(coverage_table)


#2️⃣ VARIANCE CHECK (T vs Empirical)
var_df <- plot_df %>%
  dplyr::select(Method, Parameter, EmpVar, RubinVar) %>%
  pivot_longer(cols = c(EmpVar, RubinVar),
               names_to = "Type",
               values_to = "Variance") %>%
  mutate(Type = recode(Type,
                       EmpVar = "Empirical Variance\n(Var of pooled theta across runs)",
                       RubinVar = "Rubin Total Variance\n(Average T = W + (1+1/M)B)"))

max_var <- max(var_df$Variance, na.rm = TRUE)

p_var <- ggplot(var_df,
                aes(x = Parameter,
                    y = Variance,
                    fill = Type)) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_wrap(~ Method) +
  scale_y_continuous(limits = c(0, max_var * 1.05)) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "bottom") +
  labs(title = "Empirical Variance vs Rubin Total Variance",
       x = "Parameter",
       y = "Variance",
       fill = "")

p_ratio <- ggplot(ratio_df,
                  aes(x = Parameter,
                      y = Ratio,
                      fill = Calibration)) +
  geom_col() +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             linewidth = 1) +
  facet_wrap(~ Method) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "Variance Calibration (Rubin / Empirical)",
       x = "Parameter",
       y = "T / Empirical Variance") +
  theme(legend.position = "bottom")

grid.arrange(p_var, p_ratio, ncol = 2)

variance_summary <- plot_df %>%
  group_by(Method) %>%
  summarise(
    Avg_EmpVar   = mean(EmpVar, na.rm = TRUE),
    Avg_RubinVar = mean(RubinVar, na.rm = TRUE),
    Avg_Ratio    = mean(RubinVar / EmpVar, na.rm = TRUE)
  ) %>%
  mutate(across(where(is.numeric), round, 3))

  cat("\n================ VARIANCE CALIBRATION =================\n")
  cat("Good ≈ 1 | <1 underestimation | >1 overestimation\n\n")
  
  variance_table <- plot_df %>%
    dplyr::select(Method, Parameter, EmpVar, RubinVar) %>%
    mutate(Ratio = RubinVar / EmpVar) %>%
    dplyr::select(Method, Parameter, Ratio) %>%
    pivot_wider(names_from = Method, values_from = Ratio) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  
  print(variance_table)


library(gt)

# Create gt object
variance_gt <- variance_summary %>%
  gt() %>%
  tab_header(
    title = "Variance Match",
    subtitle = "Rubin Total Variance vs Empirical Sampling Variance"
  ) %>%
  cols_label(
    Avg_EmpVar   = "Empirical Var",
    Avg_RubinVar = "Rubin Var (T)",
    Avg_Ratio    = "T / EmpVar"
  ) %>%
  fmt_number(
    columns = c(Avg_EmpVar, Avg_RubinVar, Avg_Ratio),
    decimals = 3
  ) %>%
  tab_options(
    table.font.size = px(18),
    heading.title.font.size = px(22),
    heading.subtitle.font.size = px(16)
  )

# Save vector PDF (best for LaTeX poster)
gtsave(
  variance_gt,
  "variance_match.pdf"
)

# Save high-resolution PNG
gtsave(
  variance_gt,
  "variance_match.png",
  zoom = 3
)


#3️⃣ BIAS CHECK
cat("\nBIAS CHECK\n")
cat("Good: Bias ≈ 0 | Large systematic shift = model misspecification\n\n")

p_bias <- ggplot(plot_df,
                 aes(x = Parameter,
                     y = Bias,
                     fill = Method)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, size = 26, face = "bold"),
    axis.text.y  = element_text(size = 26),
    axis.title   = element_text(size = 30),
    plot.title   = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 28),
    legend.text  = element_text(size = 26)
  ) +
  labs(title = "Bias by Parameter",
       y = "Bias",
       x = "Parameter")

print(p_bias)

ggsave("bias_plot.png",
       plot = p_bias,
       width = 12,
       height = 8,
       dpi = 600)

cat("\n================ BIAS =================\n")
cat("Good ≈ 0 | Large systematic shift = misspecification\n\n")

bias_table <- plot_df %>%
  dplyr::select(Method, Parameter, Bias) %>%
  pivot_wider(names_from = Method, values_from = Bias) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

print(bias_table)


############################################################
# RUBIN INTERNAL DIAGNOSTICS (ALL SIM RUNS)
############################################################

rubin_avg_list <- list()

m <- m_val
p <- length(true_beta)
n_obs <- nrow(full_store[[1]])
nu_com <- n_obs - p

for (method in methods) {
  
  n_runs <- length(results[[method]]$b_store)
  
  lambda_all <- matrix(NA, n_runs, p)
  fmi_all    <- matrix(NA, n_runs, p)
  df_all     <- matrix(NA, n_runs, p)
  
  for (r in seq_len(n_runs)) {
    
    B_vec  <- results[[method]]$b_store[[r]]
    T_vec  <- results[[method]]$t_store[[r]]
    SE_vec <- results[[method]]$se_store[[r]]
    
    W_vec <- SE_vec^2
    
    lambda_vec <- ((1 + 1/m) * B_vec) / T_vec
    r_vec      <- ((1 + 1/m) * B_vec) / W_vec
    
    # Barnard–Rubin df
    nu_obs <- (1 - lambda_vec) * nu_com * (nu_com + 1) / (nu_com + 3)
    nu_vec <- 1 / ((lambda_vec^2 / (m - 1)) + (1 / nu_obs))
    
    gamma_vec <- (r_vec + 2 / (nu_vec + 3)) / (1 + r_vec)
    
    lambda_all[r, ] <- lambda_vec
    fmi_all[r, ]    <- gamma_vec
    df_all[r, ]     <- nu_vec
  }
  
  B_all <- matrix(NA, n_runs, p)
  W_all <- matrix(NA, n_runs, p)
  T_all <- matrix(NA, n_runs, p)
  W_vec <- SE_vec^2
  B_all[r, ] <- B_vec
  W_all[r, ] <- W_vec
  T_all[r, ] <- T_vec
  
  rubin_avg_list[[method]] <- data.frame(
    Method = method,
    Parameter = names(true_beta),
    Mean_Lambda = colMeans(lambda_all, na.rm = TRUE),
    Mean_FMI    = colMeans(fmi_all, na.rm = TRUE),
    Mean_DF     = colMeans(df_all, na.rm = TRUE),
    Mean_B      = colMeans(B_all, na.rm = TRUE),
    Mean_W      = colMeans(W_all, na.rm = TRUE),
    Mean_T      = colMeans(T_all, na.rm = TRUE)
  )
  

rubin_avg_table <- do.call(rbind, rubin_avg_list)


rubin_method_summary <- rubin_avg_table %>%
  group_by(Method) %>%
  summarise(
    Avg_Lambda = round(mean(Mean_Lambda), 3),
    Avg_FMI    = round(mean(Mean_FMI), 3),
    Avg_DF     = round(mean(Mean_DF), 1),
    Avg_B      = round(mean(Mean_B), 4),
    Avg_W      = round(mean(Mean_W), 4),
    Avg_T      = round(mean(Mean_T), 4)
  )

}
rubin_method_summary

library(gt)

rubin_gt <- rubin_method_summary %>%
  gt() %>%
  tab_header(
    title = "Rubin Internal Diagnostics",
    subtitle = "Variance Decomposition & Properness Diagnostics"
  ) %>%
  cols_label(
    Avg_Lambda = "λ",
    Avg_FMI    = "FMI",
    Avg_DF     = "DF",
    Avg_B      = "Var (B)",
    Avg_W      = "Var (W)",
    Avg_T      = "Var (T)"
  ) %>%
  tab_options(
    table.font.size = px(18),        # good for poster
    heading.title.font.size = px(22),
    heading.subtitle.font.size = px(16)
  )
gtsave(
  rubin_gt,
  "rubin_diagnostics.png",
  zoom = 3
)
gtsave(
  rubin_gt,
  "rubin_diagnostics.pdf"
)

cat("\n================ RUBIN INTERNAL DIAGNOSTICS =================\n")
cat("Lambda = Proportion of missing-data variance\n")
cat("FMI    = Fraction of missing information\n")
cat("DF     = Barnard–Rubin degrees of freedom\n\n")

print(rubin_method_summary)
cat("==============================================================\n\n")