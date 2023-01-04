#' Estimate ATE using from AIPW estimator using cross-fit algorithm (multiple repetition)
#'
#' @param data a data frame of tibble
#' @param exposure name of exposure variable
#' @param outcome name of outcome variable
#' @param covarsT a vector of names of covaraites for treatment model
#' @param covarsO a vector of names of covaraites for outcome model
#' @param family.y it is the family for outcome model. It can `binomial() (default)` or `"gaussian"`
#' @param learners similar as \code{\link[Superlearner:SL.library()]{Superlearner::SL.library()}}
#' @param control similar as  \code{\link[Superlearner:cvControl()]{Superlearner::cvControl()}}
#' @param num_cf number of repetition done. The default is 5.
#' @param n_split number of splits used, default `n_split = 3`
#' @param rand_split logical value; if be FALSE `(default)`, discordant splits for exposure and outcome model are chosen systematically; otherwise chosen randomly.
#' @param seed numeric value to reproduce the splits distribution
#' @return a tibble of the estimates
#' @return a tibble of the estimates
#'
#' @import dplyr tibble tidyr purrr furrr tmle
#'
#' @export
#'
#' @examples
#'
#' # See the README file for details
#'
#' sum(1:4)
#'
#'
aipw_multiple_p <- function(data,
                            exposure,
                            outcome,
                            covarsT,
                            covarsO,
                            family.y = binomial(),
                            learners = c("SL.glm", "SL.glmnet", "SL.xgboost"),
                            control,
                            num_cf = 5,
                            n_split = 3,
                            rand_split = FALSE,
                            seed = 145){

  #Initialize results
  runs <- tibble(r1=double(), r0=double(), rd=double(), v1=double(), v0=double(), vd=double())

  #Run on num_cf splits
  set.seed(seed)
  cf_seed = sample(num_cf)
  for(cf in 1:num_cf){
    seed = cf_seed[cf]
    runs <- bind_rows(runs, aipw_single_p(data,
                                          exposure,
                                          outcome,
                                          covarsT,
                                          covarsO,
                                          family.y,
                                          learners,
                                          control,
                                          n_split,
                                          rand_split,
                                          seed))

  }
  #Medians of splits
  medians <- apply(runs, 2, median)


  #Corrected variance terms
  runs <- runs %>%
    mutate(mv1 = v1 + (r1-medians[1])^2,
           mv0 = v0 + (r0-medians[2])^2,
           mvd = vd + (rd-medians[3])^2)

  results <- apply(runs, 2, median)

  fit <- list()

  fit$ATE <- tibble(Estimate = results["rd"], std.error = sqrt(results["mvd"]),
                    lower_ci = results["rd"] - 1.959964*sqrt(results["mvd"]),
                    upper_ci = results["rd"] + 1.959964*sqrt(results["mvd"]))

  fit$Effct_Treat <- tibble(Estimate = results["r1"], std.error = sqrt(results["mv1"]),
                            lower_ci = results["r1"] - 1.959964*sqrt(results["mv1"]),
                            upper_ci = results["r1"] + 1.959964*sqrt(results["mv1"]))

  fit$Effct_Control <- tibble(Estimate = results["r0"], std.error = sqrt(results["mv0"]),
                              lower_ci = results["r0"] - 1.959964*sqrt(results["mv0"]),
                              upper_ci = results["r0"] + 1.959964*sqrt(results["mv0"]))

  fit

}
