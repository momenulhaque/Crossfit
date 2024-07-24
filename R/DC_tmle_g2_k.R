#'  Estimate Average Treatment Effect (ATE) using from TMLE estimator using cross-fit algorithm (generalization 2)
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
#' @param conf.level confidence limit for confidence interval, `default = 0.95`.
#' @return A tibble containing risk (or mean) difference (`ATE`), standard error (`se`), lower and upper confidence intervals (`lower.ci` and `upper.ci`, respectively).
#'
#' @section Data Generation:
#' The package includes a function, \code{gen.data}, to generate synthetic datasets for simulation studies. The function creates covariates and treatment indicators, and can generate counterfactual outcomes based on user-defined parameters.
#'
#' @references
#' Kang, J. D. Y., & Schafer, J. L. (2007). Demystifying Double Robustness: A Comparison of Alternative Strategies for Estimating a Population Mean from Incomplete Data. \emph{Statistical Science}, 22(4), 523-539.
#'
#' @import dplyr tibble tidyr purrr furrr
#'
#' @export
#'
#' @examples
#'
#' # Generate a synthetic dataset with 1200 observations
#' data <- gen.data(n = 1200, do.misp = TRUE, my.transform = TRUE, seed = 234)
#' head(data)
#'
#' # Fit the TMLE model to the generated data
#' fit_tmle_g2 <- DC_tmle_g2_k(data = data, exposure = "X", outcome = "Y",
#'                             covarsT = c("C1", "C2", "C3", "C4"),
#'                             covarsO = c("C1", "C2", "C3", "C4"),
#'                             family.y = "gaussian",
#'                             learners = c("SL.glm", "SL.mean"),
#'                             control = list(V = 3, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
#'                             num_cf = 5, n_split = 5, rand_split = FALSE,
#'                             Qbounds = 0.0025, gbounds = 0.0025, seed = 236,
#'                             conf.level = 0.95, stat = "median")
#'
#' # Display the results
#' print(fit_tmle_g2)
#'
DC_tmle_g2_k <- function(data,
                            exposure,
                            outcome,
                            covarsT,
                            covarsO,
                            family.y="binomial",
                            learners,
                            control,
                            n_split,
                            num_cf,
                            rand_split=FALSE,
                            gbounds = NULL,
                            Qbounds = 5e-04,
                            seed=146,
                            conf.level=0.95,
                            stat = "median"){

  runs <- list()
  # placeholder_output <- generate_placeholder_output(learners, n_split)
  #Run on num_cf splits
  set.seed(seed)
  cf_seed = sample(num_cf)
  n = nrow(data)
  if(is.null(gbounds)) gbounds = 5/sqrt(n)/log(n)
  #################### step 2 and 3 ########################

  for(cf in 1:num_cf){
    seed1 = cf_seed[cf]
    fit_result = try({
      tmle_single_g2_p(data,
                       exposure,
                       outcome,
                       covarsT,
                       covarsO,
                       family.y,
                       learners,
                       control,
                       n_split,
                       rand_split,
                       gbounds,
                       Qbounds,
                       seed=seed1)}, silent = TRUE)

    if (inherits(fit_result, "try-error")) {
      fit_sngle <-  data.frame(rd=NA, var = NA)
    } else {
      fit_sngle <- fit_result
    }

    runs[[cf]] <- fit_sngle
  }

  res = dplyr::bind_rows(runs)

  if(stat == "mean"){
    medians <- apply(res, 2, mean, na.rm = TRUE)
    res <- res %>% mutate(var0 = var + (rd - medians[1])^2)
    results <- apply(res, 2, mean, na.rm = TRUE)
  }
  if(stat == "median"){
    medians <- apply(res, 2, median, na.rm = TRUE)
    res <- res %>% mutate(var0 = var + (rd - medians[1])^2)
    results <- apply(res, 2, median, na.rm = TRUE)
  }
  t.value = qt((1-conf.level)/2, nrow(data), lower.tail = F)

  l_ci = results[1] - t.value*sqrt(results[3])
  u_ci = results[1] + t.value*sqrt(results[3])

  res1 = tibble(ATE=results[1], se = sqrt(results[3]), lower.ci = l_ci, upper.ci = u_ci)

  return(res1)
}
