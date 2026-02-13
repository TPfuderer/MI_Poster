library(pacman)

pacman::p_load(tidyverse, readr, dplyr, VIM, MASS, patchwork)

full_df <- read_csv("tristan stuff test/FULL_real_dataset.csv")
#View(FULL_real_dataset)
miss_df <- read_csv("tristan stuff test/analysis_dataset_with_missing.csv")
#View(analysis_dataset_with_missing)

str(miss_df)
summary(miss_df)
head(miss_df)

miss_df %>% is.na() %>% colMeans()
# or
miss_df %>% sapply(., is.na) %>% colMeans()

# left plot: Proportion of missing values
# right plot: bars shows how often variables are missing seperated or together
# blue -> observed, red -> missing 
aggr(x = miss_df,
     numbers = TRUE,          # display numbers on the right of the bars
     prop = c(TRUE, FALSE),   # left: proportion, right: total amount
     oma = c(10,5,5,3))       # outer margin area

aggr(x = miss_df, 
     combined = TRUE,         # combining the two plots into one
     numbers = TRUE,          
     prop = FALSE,
     oma = c(10,5,5,3))   

## What is the percentage of missing values in "age" for each category of "sex"?
## First variable -> splits the data -> horizontal axis of the plot
## Second variable -> interested in the proportion of missing within the subgroups
spineMiss(x = miss_df[, c("X1", "X10")])

matrixplot(miss_df, cex.axis=0.8)

par(mfrow = c(2, 5), mar = c(2,2,3,1))

matrixplot(miss_df, sortby = "X1",  cex.axis = 0.6, main = "sortby = X1")
matrixplot(miss_df, sortby = "X2",  cex.axis = 0.6, main = "sortby = X2")
matrixplot(miss_df, sortby = "X3",  cex.axis = 0.6, main = "sortby = X3")
matrixplot(miss_df, sortby = "X4",  cex.axis = 0.6, main = "sortby = X4")
matrixplot(miss_df, sortby = "X5",  cex.axis = 0.6, main = "sortby = X5")

matrixplot(miss_df, sortby = "X6",  cex.axis = 0.6, main = "sortby = X6")
matrixplot(miss_df, sortby = "X7",  cex.axis = 0.6, main = "sortby = X7")
matrixplot(miss_df, sortby = "X8",  cex.axis = 0.6, main = "sortby = X8")
matrixplot(miss_df, sortby = "X9",  cex.axis = 0.6, main = "sortby = X9")
matrixplot(miss_df, sortby = "X10", cex.axis = 0.6, main = "sortby = X10")


marginplot(miss_df[,c("X1","X2")])

