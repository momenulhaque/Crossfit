# Crossfit: An R package to fit and compare DC-AIPW and DC-TMLE for multiple splits

## How to install and apply DC-AIPW and DC-TMLE on a data set

```{r}
credentials::set_github_pat() # it asks for token. 
devtools::install_github("momenulhaque/Crossfit") # it will install the package
library(Crossfit) # Now the package is ready to use


################# Without parallelization ##################
library(tidyverse)
require(furrr)
require(tibble)
require(SuperLearner)

# I used only two learners
#Logistic Regression

SL.glm.DCDR <- function(...){
  SL.glm(...)
}

#4 degree GAM
SL.gam4.DCDR <- function(...){
  SL.gam(..., deg.gam=4)
}

# data can be from Crossfit package
df = data 

exposure="statin"
outcome="Y"


covarsT <- c("age", "ldl_log", "risk_score") # covariate for exposure
covarsO <- c("age", "ldl_log", "risk_score") # covariate for outcome

learners <- c("SL.glm.DCDR", "SL.gam4.DCDR")

control <- SuperLearner.CV.control(V=5)

## Wrapper functions

aipw_sim <- function(df, num_cf = 3, n_split, seed){
  aipw_output <- aipw_multiple_p(df, exposure, outcome, covarsT, covarsO, learners, control,
                                 num_cf, n_split, rand_split = TRUE, seed)
  return(aipw_output)
}

tmle_sim <- function(df, num_cf = 3, n_split, seed){
  tmle_output <- tmle_multiple_p(df, exposure, outcome, covarsT, covarsO, learners, control,
                                 num_cf, n_split, rand_split = TRUE, seed)
  return(tmle_output)
}



aipw_result_p3 = aipw_sim(df=data, num_cf = 3, n_split = 3, seed = 123)
tmle_result_p3 = tmle_sim(df=data, num_cf = 3, n_split = 3, seed = 123)



############### With parallelization ##############
library(parallel)
cl <- makeCluster(detectCores())

parallel::clusterEvalQ(cl, {
  
  library(Crossfit)
  library(tidyverse)
  require(furrr)
  require(tibble)
  require(SuperLearner)
  
 #Logistic Regression
  
  SL.glm.DCDR <- function(...){
    SL.glm(...)
  }
  
  #4 degree GAM
  SL.gam4.DCDR <- function(...){
    SL.gam(..., deg.gam=4)
  }
  
  df = data # data is from Crossfit package
  
  exposure="statin"
  outcome="Y"
  
  
  covarsT <- c("age", "ldl_log", "risk_score") # covariate for exposure
  covarsO <- c("age", "ldl_log", "risk_score") # covariate for outcome
  
  learners <- c("SL.glm.DCDR", "SL.gam4.DCDR")
  
  control <- SuperLearner.CV.control(V=5)
  
  
  ## Wrapper functions
  
  aipw_sim <- function(df, num_cf = 3, n_split, seed){
    aipw_output <- aipw_multiple_p(df, exposure, outcome, covarsT, covarsO, learners, control,
                                   num_cf, n_split, rand_split = TRUE, seed)
    return(aipw_output)
  }
  
  tmle_sim <- function(df, num_cf = 3, n_split, seed){
    tmle_output <- tmle_multiple_p(df, exposure, outcome, covarsT, covarsO, learners, control,
                                   num_cf, n_split, rand_split = TRUE, seed)
    return(tmle_output)
  }


})



# I just made a list of two same `data'  to apply clusterMap function.


aipw_result_p3 = clusterMap(cl, function(df, seed) aipw_sim(df, num_cf = 3, n_split = 3, seed),
                            df=list(data, data) , seed=list(1:2))


tmle_result_p3 = clusterMap(cl, function(df, seed) tmle_sim(df, num_cf = 3, n_split = 3, seed),
                            df=list(data, data) , seed=list(1:2))




```
