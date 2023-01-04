# Crossfit: An R package to apply sample splitting (cross-fit) to AIPW and TMLE in causal inference

## How to install

```{r}
devtools::install_github("momenulhaque/Crossfit") # it will install the package
library(Crossfit) 
```
Now the package is ready to use. It supports applying both AIPW and TMLE for two cases-

### Case 1: Without parallelization

 1. Install the required R packages

```{r}
library(tidyverse)
require(furrr)
require(tibble)
require(SuperLearner)
```
  2. Define the learners that you want to use in superlearner training

```{r}
#Logistic Regression
  SL.glm.DCDR <- function(...){
    SL.glm(...)
  }
  
#4 degree GAM
  SL.gam4.DCDR <- function(...){
    SL.gam(..., deg.gam=4)
  }
  
#6 degree GAM
  SL.gam6.DCDR <- function(...){
    SL.gam(..., deg.gam=6)
  }
  
#Neural Network
  SL.nnet.DCDR <- function(...){
    SL.nnet(..., size=4)
  }
  
#Random forest
  SL.randomForest.DCDR <- function(...){
    SL.randomForest(..., ntree=500, nodesize=20)
  }

#Empirical mean
  SL.mean.DCDR <- function(...){
    SL.mean(...)
  }
  

learners <- c("SL.glm.DCDR", "SL.gam4.DCDR", "SL.gam6.DCDR", "SL.nnet.DCDR", "SL.randomForest.DCDR", "SL.mean.DCDR")
```

  3. Defining the data you want to use
 
```{r}
# Read the data set that you want to use. An example data set "data" can be found in this package.
df = data 
```
  3. Defining model parameters

```{r}
exposure="statin"
outcome="Y"

covarsT <- c("age", "ldl_log", "risk_score") # covariate for exposure
covarsO <- c("age", "ldl_log", "risk_score") # covariate for outcome

# Here `V=5' indicates the number of cross-validation folds that is applied in the superlearner training.
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
```
  3. Estimating the average causal effect
  Here the parameters **num_cf** is the number of repeatation, **n_split** is the number of splits and **seed** can be usefull for comparing different methods.
```{r}
aipw_result_p3 = aipw_sim(df=data, num_cf = 3, n_split = 3, seed = 123)
tmle_result_p3 = tmle_sim(df=data, num_cf = 3, n_split = 3, seed = 123)
```

### Case 2: With parallelization
The parallelization is very usefull while simulation study is conducted for large number of times. The steps are similar to case 1, except some additional steps-

 1. Do the steps 1-3 that is described in case 1 under parallel packages. 


```{r}
############### With parallelization ##############
library(parallel)
cl <- makeCluster(detectCores())

parallel::clusterEvalQ(cl, {
  
  library(Crossfit)
  library(tidyverse)
  require(furrr)
  require(tibble)
  require(SuperLearner)
  
 # Run the codes from Step 1
 # Run the codes from Step 2
 # Run the codes from Step 3

})

```


 2. Apply the following codes for estimating average causal effect under parallel packages-
 All the data set should be stored in a list object. I just made a list of two data sets using the same data `data'  to apply clusterMap function.

```{r}
aipw_result_p3 = clusterMap(cl, function(df, seed) aipw_sim(df, num_cf = 3, n_split = 3, seed),
                            df=list(data, data) , seed=list(1:2))

tmle_result_p3 = clusterMap(cl, function(df, seed) tmle_sim(df, num_cf = 3, n_split = 3, seed),
                            df=list(data, data) , seed=list(1:2))


```



# References
Zivich PN, and Breskin A. "Machine learning for causal inference: on the use of cross-fit estimators." Epidemiology 32.3 (2021): 393-401

