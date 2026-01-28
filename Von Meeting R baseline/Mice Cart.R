#install.packages(c("pacman", "mice", "car", "tidyverse", "miceRanger", "ranger"))
# ============================================================
# MICE SIMULATION – CART IMPUTATION (same as lecture)
# ============================================================

pacman::p_load(mice, car, tidyverse)

data(SLID)
set.seed(123)

# --- specify CART imputation ---
meth_cart <- make.method(SLID)
meth_cart[c("wages", "education")] <- "cart"

# --- run MICE ---
Imp_CART <- mice(
  SLID,
  m = 10,
  maxit = 20,
  method = meth_cart,
  printFlag = FALSE
)

# --- inspect imputations ---
Imp_CART$imp %>% lapply(head)

# --- analysis ---
fit_cart <- with(
  Imp_CART,
  lm(wages ~ age + education + sex)
)

# --- pooling ---
pooled_cart <- pool(fit_cart)

# --- results ---
summary(pooled_cart)

# --- diagnostics ---
plot(Imp_CART)
densityplot(Imp_CART)

#Between imputation variance B
pooled_cart <- mice::pool(fit_cart)

B_cart <- pooled_cart$pooled$b
B_cart

