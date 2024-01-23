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
#' @param conf.level confidence limit for confidence interval, `default = 0.95`.
#' @return It return a list of two elements. The first element `ATE` is a tibble of the estimates. The `weight` is a tibble of weights of learners for two different models.
#'
#' @import dplyr tibble tidyr purrr furrr
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
DC_tmle_k <- function(data,
                      exposure,
                      outcome,
                      covarsT,
                      covarsO,
                      family.y="binomial",
                      learners,
                      control,
                      n_split,
                      num_cf,
                      rand_split,
                      gbound = 0.025,
                      alpha = 1e-17,
                      seed=146,
                      conf.level=0.95){

  runs <- list()
  placeholder_output <- generate_placeholder_output(learners, n_split)
  #Run on num_cf splits
  set.seed(seed)
  cf_seed = sample(num_cf)

  #################### step 2 and 3 ########################

  for(cf in 1:num_cf){
    seed1 = cf_seed[cf]

    fit_sngle_result <- try({
      tmle_single_p(data,
                    exposure,
                    outcome,
                    covarsT,
                    covarsO,
                    family.y,
                    learners,
                    control,
                    n_split,
                    rand_split,
                    gbound,
                    alpha,
                    seed=seed1)
    }, silent = TRUE)

    if (inherits(fit_sngle_result, "try-error")) {
      fit_sngle <- placeholder_output
    } else {
      fit_sngle <- fit_sngle_result
    }

    runs[[cf]] <- fit_sngle
  }

  res = purrr::map(runs, "results") %>%
    dplyr::bind_rows() %>%
    dplyr::filter(!is.na(rd))

  weight1 = purrr::map(runs, "weight") %>%
    dplyr::bind_rows() %>%
    dplyr::filter(!is.na(model))

  result <- weight1 %>%
    dplyr::group_by(model, split) %>%
    dplyr::summarise(across(where(is.numeric), summarize_multiple, .names = "{col}_{fn}"), .groups = "drop") %>%
    select(-matches("n.avail") | prev_n.avail)

  result_summary <- result %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(across(matches("_mean$"), mean, na.rm = TRUE, .names = "{col}"))


  medians <- apply(res, 2, median, na.rm = TRUE)

  res <- res %>% mutate(var0 = var + (rd - medians[1])^2)

  results <- apply(res, 2, median, na.rm = TRUE)

  t.value = qt((1-conf.level)/2, nrow(res), lower.tail = F)

  l_ci = results[1] - t.value*sqrt(results[3])
  u_ci = results[1] + t.value*sqrt(results[3])

  res = tibble(rd=results[1], se = sqrt(results[3]), lower.ci = l_ci, upper.ci = u_ci)

  fit <- list()

  fit$ATE = res
  fit$weight = result_summary
  # fit$weight.summary = result
  # fit$weight.all = weight1
  fit

}
