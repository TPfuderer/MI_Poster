#install.packages(c("pacman", "miceRanger", "car", "tidyverse"))
# ============================================================
# MICERANGER SIMULATION – RANDOM FOREST IMPUTATION (1:1)
# ============================================================

pacman::p_load(miceRanger, car, tidyverse)

miss_df <- read_csv("tristan stuff test/analysis_dataset_with_missing.csv")
#View(analysis_dataset_with_missing)
set.seed(123)
summary(miss_df)
# ------------------------------------------------------------
# 1) Run Random-Forest-based MICE via miceRanger
# ------------------------------------------------------------
Imp_RF <- miceRanger(
  miss_df,
  m = 25,
  maxit = 10,
  num.threads = 4,
  verbose = FALSE
)

# ------------------------------------------------------------
# 2) Inspect imputations (correct miceRanger way)
# ------------------------------------------------------------
imputed_list <- completeData(Imp_RF)

# inspect first rows of each imputed dataset
lapply(imputed_list, head)

# ------------------------------------------------------------
# 3) Analysis (equivalent to with(Imp_RF, lm(...)))
# ------------------------------------------------------------
yfit_list <- lapply(
  imputed_list,
  function(dat) lm(Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = dat)
)


# ------------------------------------------------------------
# 4) Pooling via Rubin's rules
# ------------------------------------------------------------
class(fit_list) <- "mira"
pooled_rf <- mice::pool(fit_list)

# ------------------------------------------------------------
# 5) Results (Rubin pooling output)
# ------------------------------------------------------------
summary(pooled_rf)

names(Imp_RF$data)
Imp_RF$varSummary
names(Imp_RF$imputedData)


# ------------------------------------------------------------
# 6) Diagnostics (miceRanger-native)
# ------------------------------------------------------------
imputed_vars <- rownames(Imp_RF$varSummary)
imputed_vars


# Internal miceRanger diagnostic plots
# (imputed vs observed distributions, variable-wise)
plotDistributions(
  miceObj = Imp_RF,
  vars = c("Y", paste0("X", 1:8))
)

# 1) Convergence / chain-style plot
plotVarConvergence(
  Imp_RF,
  vars = c("Y", "X1","X2","X3","X4","X5","X6","X7","X8")
)

# 2) Imputed vs observed distributions
plotDistributions(
  Imp_RF,
  vars = c("Y", "X1","X2","X3","X4","X5","X6","X7","X8")
)

# 3) Imputation uncertainty
plotImputationVariance(
  Imp_RF,
  vars = c("Y", "X1","X2","X3","X4","X5","X6","X7","X8")
)

# 4) Optional: RF fit quality
plotModelError(Imp_RF)


# ------------------------------------------------------------
# 7) Pooling + Rubin components
# ------------------------------------------------------------

class(yfit_list) <- "mira"
pooled_rf <- mice::pool(yfit_list)

summary(pooled_rf)

# Extract pooled table
pooled_table <- pooled_rf$pooled

# Between-imputation variance
B_rf <- pooled_table$b

# Within-imputation variance
W_rf <- pooled_table$ubar

# Total variance
T_rf <- pooled_table$t

B_rf
W_rf
T_rf


