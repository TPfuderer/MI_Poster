# ============================================================
# MICE SIMULATION – RANDOM FOREST IMPUTATION (CORRECT)
# ============================================================

pacman::p_load(mice, car, tidyverse)

data(SLID)
set.seed(123)

# --- specify RF imputation ---
meth_rf <- make.method(SLID)
meth_rf[c("wages", "education")] <- "rf"

# --- run MICE with Random Forest ---
Imp_RF <- mice(
  SLID,
  m = 3,
  maxit = 2,
  method = meth_rf,
  printFlag = FALSE
)

# --- inspect imputations ---
Imp_RF$imp %>% lapply(head)

# --- analysis ---
fit_rf <- with(
  Imp_RF,
  lm(wages ~ age + education + sex)
)

# --- pooling ---
pooled_rf <- pool(fit_rf)

# --- results ---
summary(pooled_rf)

# --- diagnostics ---
plot(Imp_RF)
densityplot(Imp_RF)
