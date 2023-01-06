#'  Estimate Average Treatment Effect (ATE) using from TMLE estimator using cross-fit algorithm with parallelizing
#'
#' @param data a data frame of tibble
#' @param exposure name of exposure variable
#' @param outcome name of outcome variable
#' @param covarsT a vector of names of covaraites for treatment model
#' @param covarsO a vector of names of covaraites for outcome model
#' @param family.y it is the family for outcome model. It can `binomial() (default)` or `"gaussian"`
#' @param learners similar as\code{SL.library()} in `SuperLearner` package.
#' @param control similar as \code{cvControl()} in `SuperLearner` package.
#' @param num_cf number of repetition done. The default is 5.
#' @param n_split number of splits used, default `n_split = 3`
#' @param rand_split logical value; if be FALSE `(default)`, discordant splits for exposure and outcome model are chosen systematically; otherwise chosen randomly.
#' @param gbound value between (0,1) for truncation of predicted probabilities. The defaults are 0.025 and 0.975. See \code{tmle::tmle()} for more information.
#' @param alpha used to keep predicted initial values bounded away from (0,1) for logistic fluctuation. The defaults are 1e-17 and 1-1e-17.
#' @param seed numeric value to reproduce the splits distribution
#' @return a tibble of the estimates
#'
#' @import dplyr tibble tidyr purrr furrr parallel
#'
#' @export
#'
#'
#' @examples
#'
#' # See the README file for details
#'
#' sum(1:5)
#'
#'
par_tmle <- function(data,
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
                      gbound = 0.025,
                      alpha = 1e-17,
                      seed = 145){


  cl <- makeCluster(detectCores())

  parallel::clusterExport(cl=cl,
                          varlist=c("data",
                                    "exposure",
                                    "outcome",
                                    "covarsT",
                                    "covarsO",
                                    "family.y",
                                    "learners",
                                    "control",
                                    "n_split",
                                    "rand_split",
                                    "gbound",
                                    "alpha",
                                    "seed"),
                          envir=environment()
  )


  parallel::clusterEvalQ(cl, {

  tmle_par <- function(seed, ...){

     tmle_single_p(data=data,
                   exposure=exposure,
                   outcome=outcome,
                   covarsT=covarsT,
                   covarsO=covarsO,
                   family.y=family.y,
                   learners=learners,
                   control=control,
                   n_split=n_split,
                   rand_split=rand_split,
                   gbound=gbound,
                   alpha=alpha,
                   seed)

  }

  })

  #Initialize results
  runs <- tibble(r1=double(), r0=double(), rd=double(), v1=double(), v0=double(), vd=double())

  #Run on num_cf splits
  set.seed(seed)
  cf_seed = sample(num_cf)

  runs = parallel::parLapply(cl, cf_seed, function(seed) tmle_par(seed))

  stopCluster(cl)

  runs <- dplyr::bind_rows(runs)


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

  fit$r1 <- tibble(Estimate = results["r1"], std.error = sqrt(results["mv1"]),
                   lower_ci = results["r1"] - 1.959964*sqrt(results["mv1"]),
                   upper_ci = results["r1"] + 1.959964*sqrt(results["mv1"]))

  fit$r0 <- tibble(Estimate = results["r0"], std.error = sqrt(results["mv0"]),
                   lower_ci = results["r0"] - 1.959964*sqrt(results["mv0"]),
                   upper_ci = results["r0"] + 1.959964*sqrt(results["mv0"]))

  fit

}
