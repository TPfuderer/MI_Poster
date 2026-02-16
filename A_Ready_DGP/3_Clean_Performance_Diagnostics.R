
#####################################################
# DGP END 
#########################################################
library(mice)
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
# MISLEADING PERFORMANCE METRICS (ALL SIM RUNS)
############################################################

ml_summary <- data.frame(
  Method = methods,
  Mean_RMSE_missing = NA,
  Mean_R2_missing = NA,
  Mean_SE = NA
)

for (i in seq_along(methods)) {
  
  method <- methods[i]
  
  rmse_vec <- c()
  r2_vec  <- c()
  
  for (r in seq_along(full_store)) {
    
    Imp_obj <- results[[method]]$Imp_store[[r]]
    
    if (method == "rf_ranger") {
      imputed_list <- completeData(Imp_obj)
    } else {
      imputed_list <- complete(Imp_obj, "all")
    }
    
    full_df <- full_store[[r]]
    miss_df <- miss_store[[r]]
    
    missing_index <- which(is.na(miss_df$Y))
    
    imp_y_missing <- rowMeans(
      sapply(imputed_list, function(dat) dat$Y[missing_index])
    )
    
    true_missing <- full_df$Y[missing_index]
    
    rmse_vec <- c(rmse_vec,
                  sqrt(mean((imp_y_missing - true_missing)^2)))
    
    r2_vec <- c(r2_vec,
                cor(imp_y_missing, true_missing)^2)
  }
  
  se_mat <- do.call(rbind, results[[method]]$se_store)
  
  ml_summary$Mean_RMSE_missing[i] <- mean(rmse_vec)
  ml_summary$Mean_R2_missing[i]   <- mean(r2_vec)
  ml_summary$Mean_SE[i]           <- mean(se_mat)
}

ml_summary[,-1] <- round(ml_summary[,-1],3)


cat("\n================ MISLEADING METRICS =================\n")
print(ml_summary)
cat("=====================================================\n\n")

##########
#plot
##########
p1 <- ggplot(ml_summary,
             aes(x=Method, y=Mean_MSE_missing, fill=Method)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(title="MSE on Missing Values (Lower = Better?)")

p2 <- ggplot(ml_summary,
             aes(x=Method, y=Mean_R2_missing, fill=Method)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(title="R² on Missing Values (Higher = Better?)")

p3 <- ggplot(ml_summary,
             aes(x=Method, y=Mean_SE, fill=Method)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(title="Mean Standard Error (Smaller = Better?)")

gridExtra::grid.arrange(p1, p2, p3, ncol=1)

