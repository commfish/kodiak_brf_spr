library(TMB)
library(tidyverse)
library(here)
library(funcr) #devtools::install_github("ben-williams/funcr")
theme_set(theme_report())


# scale data to 1
range01 <- function(x){
  (x - min(x)) / (max(x) - min(x))
}