# Crossfit: An R Package Apply Double Cross-fit Approach to TMLE in Causal Inference
---
Author: Momenul Haque Mondol & Mohammad Ehsanul Karim
---

## How to install

```{r}
remotes::install_github("momenulhaque/Crossfit") # it will install the package
library(Crossfit) 
```
Now the package is ready to use. It supports applying both AIPW and TMLE. 

 1. Install the required R packages

```{r}
require(SuperLearner)
```
 
 2. Defining the data you want to use
 
```{r}
# Read the data set that you want to use. An example data set "statin_sim_data" can be found in this package.
data = statin_sim_data 
```
 3. Defining the model parameters

```{r}
exposure="statin"
outcome="Y"
covarsT = c("age", "ldl_log", "risk_score") # covariate for exposure model
covarsO = c("age", "ldl_log", "risk_score") # covariate for outcome model
family.y = "binomial"
learners=c("SL.glm", "SL.glmnet", "SL.xgboost")
control=list(V = 3, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)
num_cf = 5 # number of repetitions
n_split = 4 # number of splits
rand_split = FALSE # splits' crossing pattern is not random
gbound = 0.025
alpha = 1e-17
seed = 156
conf.level = 0.95 # confidence level for confidence interval (default 0.95)
```

 4. Estimating the average treatment effect (ATE) using generalization 1. 

```{r}
fit_tmle_g1 <- DC_tmle_g1_k(data,
                        exposure,
                        outcome,
                        covarsT,
                        covarsO,
                        family.y,
                        learners,
                        control,
                        num_cf, 
                        n_split ,
                        rand_split,
                        gbound,
                        alpha,
                        seed,
                        conf.level)

```

 5. Understanding the results

The object `fit_tmle_g1` contains risk difference (`ATE`), standard error (`se`), lower and upper confidence interval (`lower.ci` and `upper.ci` respectively). 

```{r}
fit_tmle_g1

# A tibble: 1 × 4
#      ATE     se lower.ci upper.ci
#    <dbl>  <dbl>    <dbl>    <dbl>
#   -0.115 0.0151   -0.157  -0.0726

```



 6. Estimating the ATE using generalization 2. 

```{r}
fit_tmle_g2 <- DC_tmle_g2_k(data,
                        exposure,
                        outcome,
                        covarsT,
                        covarsO,
                        family.y,
                        learners,
                        control,
                        num_cf, 
                        n_split ,
                        rand_split,
                        gbound,
                        alpha,
                        seed,
                        conf.level)

```

5. Understanding the results for generalization 2

The object `fit_tmle_g2` contains risk difference (`ATE`), standard error (`se`), lower and upper confidence interval (`lower.ci` and `upper.ci` respectively). 

```{r}
fit_tmle_g2

# A tibble: 1 × 4
#      ATE     se lower.ci upper.ci
#    <dbl>  <dbl>    <dbl>    <dbl>
#   -0.114 0.0208   -0.172  -0.0566

```


# References
Zivich PN, and Breskin A. "Machine learning for causal inference: on the use of cross-fit estimators." Epidemiology 32.3 (2021): 393-401
