
#####################################################
# DGP END 
#########################################################
library(miceRanger)

#############
#Load
burgette_results <- readRDS("C:/Users/pfudi/PycharmProjects/MI_Poster/tristan stuff test/RDS_n=100.rds")

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
# RUBIN PROPERNESS DIAGNOSTICS (CLEAN VERSION)
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

# ----------------------------------------------------------
# Build Properness Table
# ----------------------------------------------------------

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
  df_mat   <- do.call(rbind, results[[method]]$df_store)
  b_mat    <- do.call(rbind, results[[method]]$b_store)
  t_mat    <- do.call(rbind, results[[method]]$t_store)
  
  R_sim <- nrow(coef_mat)
  
  colnames(coef_mat) <- names(true_beta)
  
  # ----- Bias -----
  bias_vec <- colMeans(coef_mat) - true_beta
  mean_bias <- mean(abs(bias_vec))
  
  # ----- Coverage -----
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
  
  # ----- Variance comparison -----
  emp_var <- apply(coef_mat, 2, var)
  mean_T  <- colMeans(t_mat)
  T_ratio <- mean(mean_T / emp_var, na.rm=TRUE)
  
  # Store table values
  properness_table$Mean_Bias[i] <- mean_bias
  properness_table$Mean_Coverage[i] <- mean(coverage)
  properness_table$Mean_T_over_EmpVar[i] <- T_ratio
  
  # Store plot data
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

cat("\n================= PROPERNESS SUMMARY =================\n")
print(properness_table)
cat("======================================================\n\n")

# ----------------------------------------------------------
# Combine Plot Data
# ----------------------------------------------------------

plot_df <- bind_rows(plot_storage)
ratio_df <- plot_df %>%
  mutate(Ratio = RubinVar / EmpVar,
         Calibration = case_when(
           Ratio < 0.95 ~ "Underestimation",
           Ratio > 1.05 ~ "Overestimation",
           TRUE ~ "Well calibrated"
         ))

# ==========================================================
# 1️⃣ COVERAGE CHECK
# ==========================================================

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
  theme_minimal() +
  labs(title="Coverage by Parameter",
       y="Coverage")

print(p_cov)

# ==========================================================
# 2️⃣ VARIANCE CHECK (T / Empirical)
# ==========================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

# Rename for clarity
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
  facet_wrap(~ Method) +                     # ← fixed scale
  scale_y_continuous(limits = c(0, max_var * 1.05)) +  # ← same limits
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


# ==========================================================
# 3️⃣ BIAS CHECK
# ==========================================================

cat("\nBIAS CHECK\n")
cat("Good: Bias ≈ 0 | Large systematic shift = model misspecification\n\n")

p_bias <- ggplot(plot_df,
                 aes(x=Parameter,
                     y=Bias,
                     fill=Method)) +
  geom_bar(stat="identity",
           position=position_dodge()) +
  geom_hline(yintercept=0,
             linetype="dashed") +
  theme_minimal() +
  labs(title="Bias by Parameter")

print(p_bias)

############################################################
# RUBIN DIAGNOSTICS TABLE (FIRST SIM RUN)
############################################################

r <- 1
m <- m_val

rubin_table_list <- list()

for (method in methods) {
  
  coef_r  <- results[[method]]$coef_store[[r]]
  B_vec   <- results[[method]]$b_store[[r]]
  T_vec   <- results[[method]]$t_store[[r]]
  SE_vec  <- results[[method]]$se_store[[r]]
  
  W_vec <- SE_vec^2  # within-imputation variance
  
  # --------------------------------------------------------
  # Rubin quantities
  # --------------------------------------------------------
  
  lambda_vec <- ((1 + 1/m) * B_vec) / T_vec
  r_vec      <- ((1 + 1/m) * B_vec) / W_vec
  
  # Complete-data df (approximation: n - p)
  n_obs <- nrow(full_store[[r]])
  p     <- length(true_beta)
  nu_com <- n_obs - p
  
  # Observed-data df
  nu_obs <- (1 - lambda_vec) * nu_com * (nu_com + 1) / (nu_com + 3)
  
  # Barnard–Rubin df
  nu_vec <- 1 / ( (lambda_vec^2 / (m - 1)) + (1 / nu_obs) )
  
  # Fraction of Missing Information
  gamma_vec <- (r_vec + 2 / (nu_vec + 3)) / (1 + r_vec)
  
  rubin_table_list[[method]] <- data.frame(
    Method = method,
    Parameter = names(true_beta),
    W = W_vec,
    B = B_vec,
    T = T_vec,
    Lambda = lambda_vec,
    RIV = r_vec,
    DF = nu_vec,
    FMI = gamma_vec
  )
}

rubin_table <- do.call(rbind, rubin_table_list)

rubin_table[ , 3:ncol(rubin_table)] <- round(
  rubin_table[ , 3:ncol(rubin_table)], 4
)

cat("\n================ RUBIN DIAGNOSTICS (SIM 1) =================\n")
print(rubin_table)
cat("============================================================\n\n")

############################################################
# RUBIN DIAGNOSTICS — AVERAGED OVER ALL SIM RUNS
############################################################

rubin_avg_list <- list()

for (method in methods) {
  
  lambda_all <- list()
  riv_all    <- list()
  fmi_all    <- list()
  df_all     <- list()
  
  for (r in seq_along(results[[method]]$b_store)) {
    
    m <- m_val
    
    B_vec  <- results[[method]]$b_store[[r]]
    T_vec  <- results[[method]]$t_store[[r]]
    SE_vec <- results[[method]]$se_store[[r]]
    
    W_vec <- SE_vec^2
    
    lambda_vec <- ((1 + 1/m) * B_vec) / T_vec
    r_vec      <- ((1 + 1/m) * B_vec) / W_vec
    
    # Complete-data df approx
    n_obs <- nrow(full_store[[r]])
    p     <- length(true_beta)
    nu_com <- n_obs - p
    
    nu_obs <- (1 - lambda_vec) * nu_com * (nu_com + 1) / (nu_com + 3)
    nu_vec <- 1 / ((lambda_vec^2 / (m - 1)) + (1 / nu_obs))
    
    gamma_vec <- (r_vec + 2 / (nu_vec + 3)) / (1 + r_vec)
    
    lambda_all[[r]] <- lambda_vec
    riv_all[[r]]    <- r_vec
    fmi_all[[r]]    <- gamma_vec
    df_all[[r]]     <- nu_vec
  }
  
  lambda_mat <- do.call(rbind, lambda_all)
  riv_mat    <- do.call(rbind, riv_all)
  fmi_mat    <- do.call(rbind, fmi_all)
  df_mat     <- do.call(rbind, df_all)
  
  rubin_avg_list[[method]] <- data.frame(
    Method = method,
    Parameter = names(true_beta),
    Mean_Lambda = colMeans(lambda_mat),
    Mean_RIV = colMeans(riv_mat),
    Mean_DF = colMeans(df_mat),
    Mean_FMI = colMeans(fmi_mat)
  )
}

rubin_avg_table <- do.call(rbind, rubin_avg_list)

rubin_avg_table[ , 3:ncol(rubin_avg_table)] <- round(
  rubin_avg_table[ , 3:ncol(rubin_avg_table)], 3
)

cat("\n================ RUBIN DIAGNOSTICS (ALL SIM RUNS) =================\n")
print(rubin_avg_table)
cat("====================================================================\n\n")

############################################################
# RUBIN DIAGNOSTICS — AVERAGED OVER ALL SIM RUNS (FIXED)
############################################################

rubin_avg_list <- list()

m <- m_val
p <- length(true_beta)

# use n from first dataset (constant across sims)
n_obs <- nrow(full_store[[1]])
nu_com <- n_obs - p

for (method in methods) {
  
  lambda_all <- list()
  riv_all    <- list()
  fmi_all    <- list()
  df_all     <- list()
  
  n_runs <- length(results[[method]]$b_store)
  
  for (r in seq_len(n_runs)) {
    
    B_vec  <- results[[method]]$b_store[[r]]
    T_vec  <- results[[method]]$t_store[[r]]
    SE_vec <- results[[method]]$se_store[[r]]
    
    W_vec <- SE_vec^2
    
    lambda_vec <- ((1 + 1/m) * B_vec) / T_vec
    r_vec      <- ((1 + 1/m) * B_vec) / W_vec
    
    # Barnard–Rubin components
    nu_obs <- (1 - lambda_vec) * nu_com * (nu_com + 1) / (nu_com + 3)
    nu_vec <- 1 / ((lambda_vec^2 / (m - 1)) + (1 / nu_obs))
    
    gamma_vec <- (r_vec + 2 / (nu_vec + 3)) / (1 + r_vec)
    
    lambda_all[[r]] <- lambda_vec
    riv_all[[r]]    <- r_vec
    fmi_all[[r]]    <- gamma_vec
    df_all[[r]]     <- nu_vec
  }
  
  lambda_mat <- do.call(rbind, lambda_all)
  riv_mat    <- do.call(rbind, riv_all)
  fmi_mat    <- do.call(rbind, fmi_all)
  df_mat     <- do.call(rbind, df_all)
  
  rubin_avg_list[[method]] <- data.frame(
    Method = method,
    Parameter = names(true_beta),
    Mean_Lambda = colMeans(lambda_mat, na.rm = TRUE),
    Mean_RIV = colMeans(riv_mat, na.rm = TRUE),
    Mean_DF = colMeans(df_mat, na.rm = TRUE),
    Mean_FMI = colMeans(fmi_mat, na.rm = TRUE)
  )
}

rubin_avg_table <- do.call(rbind, rubin_avg_list)

rubin_avg_table[ , 3:ncol(rubin_avg_table)] <- round(
  rubin_avg_table[ , 3:ncol(rubin_avg_table)], 3
)

cat("\n================ RUBIN DIAGNOSTICS (ALL SIM RUNS) =================\n")
print(rubin_avg_table)
cat("====================================================================\n\n")

library(dplyr)

rubin_method_summary <- rubin_avg_table %>%
  group_by(Method) %>%
  dplyr::summarise(
    Avg_Lambda = mean(Mean_Lambda),
    Avg_FMI = mean(Mean_FMI),
    Avg_DF = mean(Mean_DF)
  )

print(rubin_method_summary)

