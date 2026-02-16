
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
# STANDARD DIAGNOSTIC PLOTS — FIRST SIM RUN (WITH HEADERS)
############################################################

library(gridExtra)
library(grid)

r <- 1  # chosen simulation

############################################################
# 1️⃣ miceRanger DIAGNOSTICS
############################################################

Imp_rf <- results[["rf_ranger"]]$Imp_store[[r]]

# Convergence
p_rf_conv <- plotVarConvergence(Imp_rf) +
  labs(title = "miceRanger — Convergence") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

print(p_rf_conv)

# Distribution
p_rf_dens <- plotDistributions(Imp_rf) +
  labs(title = "miceRanger — Distribution") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

print(p_rf_dens)


############################################################
# 2️⃣ mice PMM DIAGNOSTICS
############################################################

Imp_pmm <- results[["pmm"]]$Imp_store[[r]]

# Convergence trace
plot(Imp_pmm,
     layout = c(1,1),
     main = "PMM (mice) — Convergence",
     cex.main = 1.4)

# Density plot
densityplot(Imp_pmm,
            main = "PMM (mice) — Distribution",
            cex.main = 1.4)


############################################################
# 3️⃣ mice CART DIAGNOSTICS
############################################################

Imp_cart <- results[["cart"]]$Imp_store[[r]]

# Convergence trace
plot(Imp_cart,
     layout = c(1,1),
     main = "CART (mice) — Convergence",
     cex.main = 1.4)

# Density plot
densityplot(Imp_cart,
            main = "CART (mice) — Distribution",
            cex.main = 1.4)


############################################################
# CHECK MISSINGNESS MECHANISM — FIRST SIM RUN
############################################################

r <- 1

cat("\n================ MISSINGNESS CHECK (SIM 1) =================\n")

full_df <- full_store[[r]]
miss_df <- miss_store[[r]]

# ----------------------------------------------------------
# 1️⃣ Missingness rate
# ----------------------------------------------------------

missing_rate <- mean(is.na(miss_df$Y))
cat("Missingness rate (Y):", round(missing_rate, 3), "\n\n")


# ----------------------------------------------------------
# Create missingness indicators for all variables
# ----------------------------------------------------------

vars_check <- colnames(full_df)

miss_indicator_df <- miss_df %>%
  mutate(across(all_of(vars_check),
                ~ as.integer(is.na(.)),
                .names = "miss_{.col}"))

# ----------------------------------------------------------
# Combine with drivers X9 and X10
# ----------------------------------------------------------

plot_df <- bind_cols(
  full_df[, c("X9", "X10")],
  miss_indicator_df %>% dplyr::select(starts_with("miss_"))
)

# ----------------------------------------------------------
# Reshape to long format
# ----------------------------------------------------------

plot_long <- plot_df %>%
  pivot_longer(cols = starts_with("miss_"),
               names_to = "Variable",
               values_to = "Missing") %>%
  mutate(Variable = gsub("miss_", "", Variable))

# ----------------------------------------------------------
# Bin drivers and compute missingness rate
# ----------------------------------------------------------

bins <- 20

plot_binned <- plot_long %>%
  pivot_longer(cols = c(X9, X10),
               names_to = "Driver",
               values_to = "DriverValue") %>%
  group_by(Variable, Driver,
           bin = ntile(DriverValue, bins)) %>%
  summarise(
    DriverMean = mean(DriverValue),
    MissingRate = mean(Missing) * 100,
    .groups = "drop"
  )

# ----------------------------------------------------------
# Plot
# ----------------------------------------------------------

ggplot(plot_binned,
       aes(x = DriverMean,
           y = MissingRate,
           color = Driver)) +
  geom_point(alpha = 0.7) +
  geom_smooth(se = FALSE) +
  facet_wrap(~ Variable, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Missingness per Variable vs Drivers",
    x = "Driver value",
    y = "Probability of missing (%)"
  )




