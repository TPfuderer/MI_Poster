
#####################################################
# DGP END 
#########################################################

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
