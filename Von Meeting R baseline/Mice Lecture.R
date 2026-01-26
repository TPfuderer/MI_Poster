# Test 
#install.packages(c("pacman", "mice", "car", "tidyverse"))

# ============================================================
# MICE SIMULATION – FULL WORKFLOW (1:1 from lecture, runnable)
# ============================================================

# Load required packages
pacman::p_load(mice, car, tidyverse)

# Load example data (Survey of Labour and Income Dynamics)
data(SLID)

# Set seed for reproducibility
set.seed(123)

# ------------------------------------------------------------
# 1) Inspect raw data
# ------------------------------------------------------------
head(SLID)
summary(SLID)

# ------------------------------------------------------------
# 2) Run MICE (vanilla setup, reduced m and iterations)
# ------------------------------------------------------------
ImpMice <- mice(
  SLID,
  m = 15,        # number of imputations
  maxit = 30     # number of MCMC iterations
)

# ------------------------------------------------------------
# 3) Inspect imputed values
# ------------------------------------------------------------
# Show first imputed values per variable
ImpMice$imp %>% lapply(head)

# ------------------------------------------------------------
# 4) Fit model on each imputed dataset
# ------------------------------------------------------------
fit <- with(
  ImpMice,
  lm(wages ~ age + education + sex)
)

# ------------------------------------------------------------
# 5) Pool results (Rubin's rules)
# ------------------------------------------------------------
pooledRes <- pool(fit)

# ------------------------------------------------------------
# 6) Inspect pooled estimates
# ------------------------------------------------------------
summary(pooledRes)

# ------------------------------------------------------------
# 7) Diagnostics (optional but shown in lecture)
# ------------------------------------------------------------
plot(ImpMice)          # convergence plots
densityplot(ImpMice)  # observed vs imputed distributions
