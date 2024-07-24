# Crossfit: An R Package to Apply Double Cross-fit Approach to TMLE in Causal Inference
---
Author: Momenul Haque Mondol & Mohammad Ehsanul Karim
---
## How to install

```{r}
remotes::install_github("momenulhaque/Crossfit") # it will install the package
library(Crossfit) 
```
Now the package is ready to use. 

## Install other required R packages
The package primarily requires 'SuperLearner'. Additionally, it requires all the R packages associated with advanced machine learning algorithms.
```{r}
require(SuperLearner)
# require(randomForest) for SL.randoForest.
```
## Example data set 
We generated an example data set ("ObsData") using Kang at al. (2007) and attached it to the package.

```{r}

## Defining Learners
SL.xgboost.dcTMLE <- function(...){
  SL.xgboost(..., ntrees = 500, max_depth = 4, shrinkage = 0.1, minobspernode = 30)}
SL.glm.dcTMLE <- function(...){
  SL.glm(...)}
SL.gam4.dcTMLE <- function(...){
  SL.gam(..., deg.gam=4)}
SL.mean.dcTMLE <- function(...){
  SL.mean(...)}
``` 
### Estimating the average treatment effect (ATE) using generalization 1
```{r}
fit_tmle_g1 = DC_tmle_g1_k(data=ObsData,
                           exposure="X",
                           outcome="Y",
                           covarsT=c("C1", "C2", "C3", "C4"),
                           covarsO=c("C1", "C2", "C3", "C4"),
                           family.y="gaussian",
                           learners=c("SL.glm.dcTMLE", "SL.mean.dcTMLE", "SL.gam4.dcTMLE", "SL.xgboost.dcTMLE"),
                           control=list(V = 3, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                           num_cf=5,
                           n_split=3,
                           rand_split=FALSE,
                           Qbounds = 0.0025,
                           gbounds = 0.0025,
                           seed=236,
                           conf.level=0.95,
                           stat = "median")
> fit_tmle_g1
# A tibble: 1 × 4
    ATE    se lower.ci upper.ci
  <dbl> <dbl>    <dbl>    <dbl>
1  6.61 0.491     5.64     7.57
```
The object `fit_tmle_g1` contains risk (or mean) difference (`ATE`), standard error (`se`), lower and upper confidence interval (`lower.ci` and `upper.ci` respectively). 

### Estimating the average treatment effect (ATE) using generalization 2

```{r}

fit_tmle_g2 = DC_tmle_g2_k(data=ObsData,
                           exposure="X",
                           outcome="Y",
                           covarsT=c("C1", "C2", "C3", "C4"),
                           covarsO=c("C1", "C2", "C3", "C4"),
                           family.y="gaussian",
                           learners=c("SL.glm.dcTMLE", "SL.mean.dcTMLE", "SL.gam4.dcTMLE", "SL.xgboost.dcTMLE"),
                           control=list(V = 3, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                           num_cf=5,
                           n_split=5,
                           rand_split=FALSE,
                           Qbounds = 0.0025,
                           gbounds = 0.0025,
                           seed=236,
                           conf.level=0.95,
                           stat = "median")
                           
> fit_tmle_g2
# A tibble: 1 × 4
    ATE    se lower.ci upper.ci
  <dbl> <dbl>    <dbl>    <dbl>
1  6.62 0.519     5.60     7.63                           
             

```

The object `fit_tmle_g2` contains risk (or mean) difference (`ATE`), standard error (`se`), lower and upper confidence interval (`lower.ci` and `upper.ci` respectively). 

# References
Kang, J. D. Y., & Schafer, J. L. (2007). Demystifying Double Robustness: A Comparison of Alternative Strategies for Estimating a Population Mean from Incomplete Data. Statistical Science, 22(4). https://doi.org/10.1214/07-STS227

Karim M, Mondol M. "Finding the Optimal Number of Splits and Repetitions in Double Cross-fitting Targeted Maximum Likelihood Estimators." (Under review)

Mondol M, Karim M. "Towards Robust Causal Inference in Epidemiological Research: Employing Double Cross-fit TMLE in Right Heart Catheterization Data." (Under review)


# Funding
Natural Sciences and Engineering Research Council of Canada (NSERC) Discovery Grant (PG#: 20R01603, PI: Karim) as well as UBC Work Learn program (Summer 2023, PI: Karim).

# Acknowledgements
This research was partially supported by the computational resources and services provided by Advanced Research Computing at the University of British Columbia. We extend our sincere gratitude to Paul N. Zivich and Alexander Breskin for their R codes available at [this GitHub repository](https://github.com/pzivich/publications-code/tree/master/DoubleCrossFit).

