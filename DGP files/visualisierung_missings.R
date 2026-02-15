set.seed(123)

library(MASS)  
library(mice)  
library(pacman)
pacman::p_load(plot3D,dplyr,rgl,tidyr)


# Parameter 1:1 Burgette --------------------------------------------------

n  <- 5000 # Anzahl der Beobachtungen
p  <- 10 # Anzahl der Features (X1-X10)

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

# Hilfsfunktionen ---------------------------------------------------------
# rang der Matrix goes is, Kovarianzmatrix für die Datengenerierung goes out
generate_cor <- function(p) {
  Sigma <- diag(p)
  Sigma[1:4, 1:4] <- 0.5
  diag(Sigma[1:4, 1:4]) <- 1
  # Eine Korrelation unter 1 garantiert, dass wir einen MAR-Prozess haben
  Sigma[5:10, 5:10] <- 0.3 
  diag(Sigma[5:10, 5:10]) <- 1

  Sigma
}

# Kovarianzmatrix von generate_cor goes in. Der Datensatz dem die Missings hinzugefügt werden
# goes out
generate_data <- function(n, p, Sigma, beta) {
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
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

  data.frame(Y, X)
}

# Wahrscheinlichkeit für die Missing.
logit_p <- function(z) 1 / (1 + exp(-z))
# Das argument wird z wird dann von X9 und X10 abhängen
########################################

# Dataframe von generate_data goes in. DF mit NAs gehen raus
add_missingness <- function(full_df, n) {
  # Höheres X9 => Mehr NA
  # Höheres X10 => Weniger NA
  # Spiegelt sich in der Heatmap wieder
  p_miss <- logit_p(0.5 * full_df$X9 - 0.5 * full_df$X10)

  analysis_df <- full_df
  # Hier wird gesteuert, welche Zeilen NAs haben sollen (Y und X1-X8)
  vars_with_missing <- c("Y", paste0("X", 1:8))
  # Jede Zeile wird durchlaufen und einen Vektor aus, der als
  # Maske dient. Darin enthalten 0 und 1 Einträge abhängig von einer
  # Binomverteilung. Bei 1 wird der korrespondierende Eintrag in der Dataframezeile zu einem NA
  for (v in vars_with_missing) {
    miss_indicator <- rbinom(n, 1, p_miss)
    analysis_df[miss_indicator == 1, v] <- NA
  }
  analysis_df
}

# Daten generieren und Missingness hinzufügen
Sigma   <- generate_cor(p)
full_df <- generate_data(n, p, Sigma, beta)
df      <- add_missingness(full_df, n)



bin_completeness <- function(df, x9_breaks, x10_breaks,
                             cols = c(paste0("X", 1:8), "Y")) {

  df2 <- df %>%
    mutate(
      # erstellt zwei neue Spalten x9_bin und x10_bin in denen ein Intervall gespeichert ist,
      # in dem der X9- bzw. X10- Wert ein element von ist
      # Wenn durch breaks z.b. Intervallgrenzen von [0.671,0.881) entstanden sind und 
      # X9 den Wert 0.86 hat, wird x9_bin genau dieses Intervall zugeordnet.
      x9_bin  = cut(X9,  breaks = x9_breaks,  include.lowest = TRUE, right = FALSE), 
      x10_bin = cut(X10, breaks = x10_breaks, include.lowest = TRUE, right = FALSE)
    ) %>%
    # pro Zeile: wie viele der 9 Felder haben kein NA?? Wenn keine NAs da sind natürlich 9
    rowwise() %>%
    mutate(exist_count = sum(!is.na(c_across(all_of(cols))))) %>%
    # ungroup scheint bei einer rowwise funkion notwendig zu sein
    ungroup()

  agg <- df2 %>%
    #Das filtern ist eigentlich bei unserer Datenstruktur nicht notwendig
    filter(!is.na(x9_bin), !is.na(x10_bin)) %>%
    #group_by verändert vorerst nichts wenn man das printed, aber sie sind schon logisch gruppiert
    group_by(x9_bin, x10_bin) %>%
    # Summarise gibt alle Kombinationen an quadratischen bins aus in denen mindestens ein Wert liegt.
    # Den Bins zugeordneten Werten (z.b. Anzahl der Werte) können als zusätzliche Befehle angegeben werden
    summarise(
      #n() gibt die Anzahl der Elemente innerhalb einer Gruppe aus
      n = n(),
      # sum_exist ist die Anzahl aller (X1-X8 und Y)-Werte die nicht NA sind
      # Die anzahl der existierenden Werte, oder der NA ist allerdings nicht
      # aussagekräftig, da wenn in einem quadrat bin aufgrund der Verteilung
      # mehr Werte sind, man auch bei MCAR mehr NAs an dieser Stelle erwarten kann
      sum_exist = sum(exist_count),
      # Daher wird es ins Verhältnis gesetzt. Bei MCAR würde man, das der Anteil
      # an NAs in jedem Quadrat-bin (von statistischen Schwankungen abgesehen) konstant bleibt
      # bei MAR ist der Anteil an NAs von Ort (X1,X2) abhängig. 
      pct_exist = 100 * sum_exist / (n * length(cols)),
      # durch "drop" wird das ganze zu einem normalen Dataframe
      .groups = "drop"
    )

  # Matrix für hist3D bauen (Rows = X10, Cols = X9)
  z_pct <- xtabs(pct_exist ~ x10_bin + x9_bin, agg) |> as.matrix()
  z_sum <- xtabs(sum_exist ~ x10_bin + x9_bin, agg) |> as.matrix()
  z_n   <- xtabs(n ~ x10_bin + x9_bin, agg) |> as.matrix()

  # Mitten der Bins (für Achsen)
  x9_mids  <- (x9_breaks[-1]  + x9_breaks[-length(x9_breaks)]) / 2
  x10_mids <- (x10_breaks[-1] + x10_breaks[-length(x10_breaks)]) / 2

  list(agg = agg, z_pct = z_pct, z_sum = z_sum, z_n = z_n, x9_mids = x9_mids, x10_mids = x10_mids, df2 = df2)
}

# Breaks: theoretischer Bereich statt datenabhängig, damit alle Simulationsläufe abgedeckt sind
x9_breaks  <- seq(-5, 5, length.out = 150)
x10_breaks <- seq(-5, 5, length.out = 150)



#Alle relevanten Variablen für einen Simulationsdurchlauf goes in und die Daten um einen Graph zu
# erzeugen goes out
one_run <- function(n, p, beta, x9_breaks, x10_breaks, cols = c(paste0("X",1:8),"Y")) {
  Sigma <- generate_cor(p)
  full_df <- generate_data(n, p, Sigma, beta)
  df <- add_missingness(full_df, n)

  res <- bin_completeness(df, x9_breaks, x10_breaks, cols = cols)

  # res$z_pct ist eine Matrix: rows = x10 bins, cols = x9 bins
  # Leere Bins (keine Beobachtungen) auf NA setzen statt 0,
  # damit sie den Durchschnitt nicht verzerren
  z <- res$z_pct
  z[res$z_n == 0] <- NA
  z
}


# one_run wird k-mal ausgeführt und die Prozente anschließend über alle k durschläufe gemittelt
simulate_mean_matrix <- function(k, n, p, beta, x9_breaks, x10_breaks) {
  z_list <- vector("list", k)

  for (i in 1:k) {
    z_list[[i]] <- one_run(n, p, beta, x9_breaks, x10_breaks)
  }

  # Mittelwert elementweise über alle Matrizen (NA-aware, damit leere Bins
  # den Durchschnitt nicht nach unten ziehen)
  z_array <- simplify2array(z_list)
  z_mean  <- apply(z_array, 1:2, mean, na.rm = TRUE)
  z_mean
}


# Ausführung des Codes
z_mean <- simulate_mean_matrix(k = 300, n = n, p = p, beta = beta,
                               x9_breaks = x9_breaks, x10_breaks = x10_breaks)


# Mittelpunkte der Rasterquadrate berechnen für den 3d Plot, oder Heatmap
x9_mids  <- (x9_breaks[-1]  + x9_breaks[-length(x9_breaks)]) / 2
x10_mids <- (x10_breaks[-1] + x10_breaks[-length(x10_breaks)]) / 2

plot3D::hist3D(x = x10_mids, y = x9_mids, z = z_mean,
  xlab = "X10", ylab = "X9", zlab = "Ø Existenz (%)",
  col = jet.col(100),
  shade = 0.3,          
  border = NA,          
  alpha = 0.8,          

  theta = 230, phi = 20)


 image2D(z = z_mean, x = x10_mids, y = x9_mids,
  xlab = "X10", ylab = "X9", main = "Ø Existenz der Einträge in (%)",
  col = hcl.colors(100, "YlOrRd", rev = TRUE),
  cex.lab  = 2, cex.axis = 1.9, cex.main = 1.8,
  colkey = list(cex.axis = 1.3, cex.clab = 1.5))