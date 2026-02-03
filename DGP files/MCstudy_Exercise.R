set.seed(123)

library(mice)

# number of simulation cycles
R <- 50

# store the results (rel. bias and coverage) in a matrix for each method
resBD <- matrix(nrow = R, ncol = 4)
colnames(resBD) <- c("relBias_mean", "Coverage_mean", 
                     "relBias_quant","Coverage_quant")

resMIpmm <- resMIcart <- resMIrf <- resBD
# -> we get for each method a matrix with R iterations as rows and for each
#    iteration we have the rel. bias and coverage for the mean of y as columns
head(resBD)

# number of multiple imputations
M <- 10

# number of iterations for the chained equations (see mice ALGO: t = 1, ...,T)
niter <- 1

# sample size
n <- 2000

## functions for diagnostics

coverage <- function(true.value, CI.low, CI.up){
  ifelse(test = CI.low <= true.value && CI.up >= true.value, yes = 1, no = 0)
}

rel.bias <- function(true.value, est) {
  rel_bias <- (est - true.value) / true.value
  return(rel_bias)
}

## function for Rubin's combining rules, gives pooled estimates, CI

MI_analysis <- function(Q.hat,U.hat,m){
  
  if (class(Q.hat)!="numeric") {
    stop("Estimator vector for all imputations has a different class")
  }
  else{
    # pooled estimator
    Q_bar <- sum(Q.hat)/m
    # within-variance
    U_bar <- sum(U.hat)/m
    # between-variance
    B <- (1/(m-1))*sum((Q.hat-Q_bar)^2)
  }
  
  # total variance
  Tot <- U_bar+B+B/m
  # degrees of freedom
  df <- (m-1)*(1+(m/(m+1))*U_bar/B)^2
  
  # confidence intervals
  CIlow <- Q_bar-qt(0.975,df)*sqrt(Tot)
  CIup <- Q_bar+qt(0.975,df)*sqrt(Tot)
  r <- (B+B/m)/U_bar
  
  # fraction of missing information
  FMI <- (r+2/(df+3))/(1+r)
  
  # t-test
  t_value <- Q_bar/sqrt(Tot)
  
  # p-value
  p.value <-2*(1-pt(abs(t_value),df))
  
  return(cbind(Q_bar,CIlow,CIup))
}

quantVar <- function(y, p = 0.9){
  q <- as.numeric(quantile(y, probs = p))
  
  dens <- density(y, na.rm = TRUE)
  f_q <- approx(dens$x, dens$y, xout = q)$y
  
  sigma <- (p * (1 - p)) / (f_q^2)
  return(sigma)
}


for(r in 1:R) {
  
  ### Data generation
  
  # one normally and one uniformly distributed variable
  x1 <- rnorm(n = n, mean = 0, sd = 1)
  x2 <- runif(n = n, min = 0, max = 2)
  
  # continuous outcome variable with non-linear terms and non-linear error term
  eps <- rchisq(n = n, df = 107/96)
  y <- -3 + 3.5 * x1^2 + 2.75 * x2^3 + (eps - 107/96)
  
  data <- as.data.frame(cbind(x1, x2, y))
  data_bd <- data
  
  # TRUE VALUES (before deletion)
  true_mean  <- mean(data_bd$y)
  true_quant <- as.numeric(quantile(data_bd$y, 0.9))
  
  ### Before deletion
  
  # mean
  mean_est <- mean(data_bd$y)
  emp_sd <- sd(data_bd$y)
  
  ci_up <- mean_est + (qnorm(p = 0.975)*emp_sd)/sqrt(length(data_bd$y))
  ci_low <- mean_est - (qnorm(p = 0.975)*emp_sd)/sqrt(length(data_bd$y))
  
  # quantile 
  quant_est <- quantile(data_bd$y, probs = 0.9)
  quant_sd <- sqrt(quantVar(y =data_bd$y, p = 0.9)) 
  
  ci_up_quant <- quant_est + (qnorm(p = 0.975)*quant_sd)/sqrt(length(data_bd$y))
  ci_low_quant <- quant_est - (qnorm(p = 0.975)*quant_sd)/sqrt(length(data_bd$y))
  
  
  resBD[r,] <- c(rel.bias(true_mean, mean_est),
                 coverage(true_mean, ci_low, ci_up),
                 rel.bias(true_quant, quant_est),
                 coverage(true_quant, ci_low_quant, ci_up_quant))
  
  
  
  ### Generation of missing data
  
  # MAR in y
  z1 <-rnorm(n = n, mean = 0, sd = sqrt(16)) * x2
  lin_y <- 1.75 - 1.5 * z1
  prob_y <- pnorm(lin_y)    # yields around 35% missing values
  res_y <- rbinom(n = n, size = 1, prob = prob_y) 
  
  data$y[which(res_y==0)] <- NA
  
  ### Multiple imputation and analysis
  
  # use all variables as predictors, but just consider linear relationships..
  
  ## initialization
  ini <- mice(data,maxit=0)
  
  ## PMM (default method)
  imp_pmm <- mice(data, m=M, maxit=niter, print=FALSE)
  
  mice_pmm <- complete(imp_pmm, action="long", include=FALSE)
  
  # mean for each imputed data set
  theta_MI <- aggregate(y ~.imp, data=mice_pmm, mean)$y
  #theta_MI <- mice_pmm %>% group_by(.imp) %>% summarise(mean(y))
  # variance of mean for each imputed data set
  var_MI <-   aggregate(y ~.imp, data=mice_pmm, var)$y/n
  
  quant_MI <- aggregate(y ~.imp, data = mice_pmm, quantile, 0.9)$y
  quant_MI_var <- aggregate(y ~.imp, data = mice_pmm, quantVar, p = 0.9)$y/n
  
  
  ana_pmm <- MI_analysis(theta_MI,var_MI,m = M)   # -> Q_bar,CIlow,CIup
  ana_pmm_quant <- MI_analysis(quant_MI, quant_MI_var,m = M)
  
  resMIpmm[r,] <- c(rel.bias(true_mean, ana_pmm[1]), 
                    coverage(true_mean, ana_pmm[2], ana_pmm[3]),
                    rel.bias(true_quant, ana_pmm_quant[1]),
                    coverage(true_quant, ana_pmm_quant[2], ana_pmm_quant[3]))
  
  
  ## cart
  meth <- ini$method
  meth["y"] <- "cart"
  
  imp_cart <- mice(data, m=M, maxit=niter, method=meth, print=FALSE)
  mice_cart <- complete(imp_cart, action="long", include=FALSE)
  
  # mean for each imputed data set
  theta_MI <- aggregate(y ~.imp, data=mice_cart, mean)$y
  # variance of mean for each imputed data set
  var_MI <-   aggregate(y ~.imp, data=mice_cart, var)$y/n
  
  quant_MI <- aggregate(y ~.imp, data = mice_cart, quantile, 0.9)$y
  quant_MI_var <- aggregate(y ~.imp, data = mice_cart, quantVar, p = 0.9)$y/n
  
  
  ana_cart <- MI_analysis(theta_MI,var_MI,M)   # -> Q_bar,CIlow,CIup
  ana_cart_quant <- MI_analysis(quant_MI, quant_MI_var,m = M)
  
  resMIcart[r,] <- c(rel.bias(true_mean, ana_cart[1]), 
                     coverage(true_mean, ana_cart[2], ana_cart[3]),
                     rel.bias(true_quant, ana_cart_quant[1]),
                     coverage(true_quant, ana_cart_quant[2], ana_cart_quant[3]))
  
  ## random forest
  meth <- ini$method
  meth["y"] <- "rf"
  
  imp_rf <- mice(data, m=M, maxit=niter, method=meth, print=FALSE)
  mice_rf <- complete(imp_rf, action="long", include=FALSE)
  
  # mean for each imputed data set
  theta_MI <- aggregate(y ~.imp, data=mice_rf, mean)$y
  # variance of mean for each imputed data set
  var_MI <-   aggregate(y ~.imp, data=mice_rf, var)$y/n
  
  quant_MI <- aggregate(y ~.imp, data = mice_rf, quantile, 0.9)$y
  quant_MI_var <- aggregate(y ~.imp, data = mice_rf, quantVar, p = 0.9)$y/n
  
  ana_rf <- MI_analysis(theta_MI,var_MI,M)   # -> Q_bar,CIlow,CIup
  ana_rf_quant <- MI_analysis(quant_MI, quant_MI_var,m = M)
  
  resMIrf[r,] <- c(rel.bias(true_mean, ana_rf[1]), 
                   coverage(true_mean, ana_rf[2], ana_rf[3]),
                   rel.bias(true_quant, ana_rf_quant[1]),
                   coverage(true_quant, ana_rf_quant[2], ana_rf_quant[3]))
}

Bias <- Coverage <- matrix(nrow=2, ncol=4)
colnames(Bias) <- colnames(Coverage) <- c("BD","pmm","cart","rf")
rownames(Bias) <- rownames(Coverage) <- c("mean(y)","quantile")

Bias[, "BD"] <- colMeans(resBD)[c(1,3)]
Bias[, "pmm"] <- colMeans(resMIpmm)[c(1,3)] 
Bias[, "cart"] <- colMeans(resMIcart)[c(1,3)] 
Bias[, "rf"] <- colMeans(resMIrf)[c(1,3)] 

Coverage[, "BD"] <- colMeans(resBD)[c(2,4)]
Coverage[, "pmm"] <- colMeans(resMIpmm)[c(2,4)] 
Coverage[, "cart"] <- colMeans(resMIcart)[c(2,4)] 
Coverage[, "rf"] <- colMeans(resMIrf)[c(2,4)] 

round(Bias, 4)

round(Coverage, 4)