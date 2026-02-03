#install.packages(c("pacman", "miceRanger", "car", "tidyverse"))
# ============================================================
# MICERANGER SIMULATION – RANDOM FOREST IMPUTATION (1:1)
# ============================================================

pacman::p_load(miceRanger, car, tidyverse)

data(SLID)
set.seed(123)

# ------------------------------------------------------------
# 1) Run Random-Forest-based MICE via miceRanger
# ------------------------------------------------------------
Imp_RF <- miceRanger(
  SLID,
  m = 10,
  maxit = 20,
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
fit_list <- lapply(
  imputed_list,
  function(dat) lm(wages ~ age + education + sex, data = dat)
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


# ------------------------------------------------------------
# 6) Diagnostics (miceRanger-native)
# ------------------------------------------------------------

# Internal miceRanger diagnostic plots
# (imputed vs observed distributions, variable-wise)
plotDistributions(
  miceObj = Imp_RF,
  vars = c("wages", "education")
)

# 1) Convergence / chain-style plot
plotVarConvergence(
  Imp_RF,
  vars = c("wages", "education")
)

# 2) Imputed vs observed distributions
plotDistributions(
  Imp_RF,
  vars = c("wages", "education")
)

# 3) Imputation uncertainty
plotImputationVariance(
  Imp_RF,
  vars = c("wages", "education")
)

# 4) Optional: RF fit quality
plotModelError(Imp_RF)


# Checking between chain
# --- miceRanger fits ---
class(fit_list) <- "mira"
pooled_rf <- mice::pool(fit_list)

# Get Rubin components (W, B, T)
summ_rf <- summary(pooled_rf)

# B: between-imputation variance per coefficient
B_rf <- pooled_rf$pooled$b
B_rf
