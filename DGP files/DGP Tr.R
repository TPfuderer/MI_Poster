# Prozess von Burgette & Reiter 2010
# Inspiration von little sim study.html (Übung)
# MAR 
# n = 1000
# Goal = properness 

set.seed(123)

library(MASS)  
library(mice)   


# Parameter 1:1 Burgette --------------------------------------------------

n  <- 1000
p  <- 10

beta <- c(
  b0 = 0,
  b1 = 0.5,
  b2 = 0.5,
  b3 = 0.5,
  b4 = 0.5,
  b5 = 0.5,
  b6 = 1,    # X3^2
  b7 = 1,    # X1*X2
  b8 = 1     # X8*X9
)


# 1.2 Korrelation zwischen Variablen --------------------------------------

Sigma <- diag(p)

# First block: X1–X4 (corr = 0.5)
Sigma[1:4, 1:4] <- 0.5
diag(Sigma[1:4, 1:4]) <- 1

# Second block: X5–X10 (corr = 0.3)
Sigma[5:10, 5:10] <- 0.3
diag(Sigma[5:10, 5:10]) <- 1


# 2. X Y generierung (Interaktionen etc ----------------------------------

X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
colnames(X) <- paste0("X", 1:10)

epsilon <- rnorm(n)

Y <- beta["b0"] +
  beta["b1"] * X[, "X1"] +
  beta["b2"] * X[, "X2"] +
  beta["b3"] * X[, "X3"] +
  beta["b4"] * X[, "X8"] +
  beta["b5"] * X[, "X9"] +
  beta["b6"] * X[, "X3"]^2 +
  beta["b7"] * X[, "X1"] * X[, "X2"] +
  beta["b8"] * X[, "X8"] * X[, "X9"] +
  epsilon

FULL_real_dataset <- data.frame(Y, X)

# 2.3 Missingness erzeugen ------------------------------------------------

logit_p <- function(z) 1 / (1 + exp(-z))

p_miss <- logit_p(0.5 * FULL_real_dataset$X9 -
                    0.5 * FULL_real_dataset$X10)

# Nur auf X1 bis X8 anwenden (X9-10 = full)

analysis_dataset <- FULL_real_dataset

vars_with_missing <- c("Y", paste0("X", 1:8))

for (v in vars_with_missing) {
  miss_indicator <- rbinom(n, 1, p_miss)
  analysis_dataset[miss_indicator == 1, v] <- NA
}

# Unter 25% Komplette Zeilen
# ca 17% Missingness pro Variable

# ===============================
# FINAL OUTPUT FOR SIMULATION
# ===============================

FULL_real_dataset      # Oracle data (true DGP, no missingness)
analysis_dataset       # MAR data used for MI

# ============================================================
# CSV EXPORT (Ground truth + MAR data)
# ============================================================

write.csv(
  FULL_real_dataset,
  file = "FULL_real_dataset.csv",
  row.names = FALSE
)

write.csv(
  analysis_dataset,
  file = "analysis_dataset_with_missing.csv",
  row.names = FALSE
)





