install.packages("pacman")
pacman::p_load("VIM","rgl","MASS","tidyverse") 

set.seed(1)

n <- 250
X1 <- rnorm(n, 8, 3)
X2 <- 10 - 0.5 * X1 + rnorm(n, 0, 3)
X3 <- 5 + 0.6 * X1 + 0.5 * X2 + rnorm(n, 0, sqrt(2))
data1 <- as.data.frame(cbind(X1, X2, X3))

# Quantities of interest for the original sample data set
(bef.Imp <- cbind(mean = mean(X3),
                  var = var(X3),
                  corX1_X3 = cor(X1, X3))
)

X3NonNA <- X3 # save for later

# Generate missing values
data1$X3[sample(x = 1:n, size = n/2, replace = FALSE)] <- NA

# Missing indicators for missing data in X3
misind <- is.na(data1$X3)

# Indicators for observed (not missing) data in X3
obsind <- !is.na(data1$X3)

layout(matrix(c(1,2,2), nrow=1))
barMiss(data1, pos=1)
marginplot(data1[,-2])

theme_set(theme_bw()) # set global theme for ggplot2
#MisInd <- 242 #chosen randomly
set.seed(42)
MisInd <- sample(which(is.na(data1$X3)),size = 3) #chose 3 random missings to display

data1Comp <- data1 %>% mutate("X3Comp"=X3NonNA)
ggplot(data = data1)+
  geom_point(aes(x = X1, y = X3))+
  geom_point(data = data1Comp[MisInd, ], 
             aes(x = X1, y = X3Comp), 
             col ="red", size = 4, shape = 8)+
  annotate(geom = "point", x = 0, y = 10, shape = 8, size = 4, col = "red") + 
  annotate(geom = "text", x = 0.2, y = 10, label = "missing observation", hjust = "left")