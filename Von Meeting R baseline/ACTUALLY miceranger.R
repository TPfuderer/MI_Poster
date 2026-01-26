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
  m = 3,
  maxit = 2,
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

