#' Generate Synthetic Data
#'
#' The `gen.data` function generates synthetic datasets for simulation studies.
#' It can produce observed data and counterfactual outcomes based on specified parameters.
#'
#' @param n Integer. The number of observations to generate.
#' @param do.pop Logical. If `TRUE`, the function returns the propensity score and counterfactual outcomes.
#' @param do.misp Logical. If `TRUE`, the function applies a transformation to the covariates.
#' @param my.transform Logical. If `TRUE` and `do.misp` is `TRUE`, the function uses the Balzer & Westling transformation. Otherwise, it uses the Naimi et al. transformation.
#' @param seed Integer. The random seed for reproducibility.
#'
#' @details
#' The function generates four covariates (C1, C2, C3, C4) from standard normal distributions.
#' A propensity score is calculated using these covariates, and a binary treatment indicator (X) is sampled.
#' The outcome (Y) is generated based on the covariates and treatment.
#'
#' If `do.pop` is `TRUE`, the function returns the propensity score and counterfactual outcomes (Y1, Y0).
#' If `do.misp` is `TRUE`, it applies transformations to the covariates before generating the outcome.
#'
#' @return A data frame with the generated data. If `do.pop` is `TRUE`, it contains the propensity score and counterfactual outcomes.
#' If `do.misp` is `TRUE`, the covariates may be transformed.
#'
#' @examples
#' # Generate a dataset with 100 observations, including counterfactual outcomes
#' data <- gen.data(n = 100, do.pop = TRUE)
#'
#' # Generate a dataset with transformed covariates
#' data_misp <- gen.data(n = 100, do.misp = TRUE, my.transform = TRUE)
#'
#' @export
gen.data <- function(n, do.pop=F, do.misp=F, my.transform=F, seed=123){
  set.seed(seed)
  C1 <- rnorm(n,0,1)
  C2 <- rnorm(n,0,1)
  C3 <- rnorm(n,0,1)
  C4 <- rnorm(n,0,1)

  pscore <- plogis(-1 + log(1.75) * (C1 + C2 + C3 + C4))
  X <- rbinom(n, 1, pscore)

  eps <- rnorm(n,0,6)
  get.Y <- function(C1, C2, C3, C4, X, eps){
    120 + 6 * X + 3 * (C1 + C2 + C3 + C4) + eps
  }
  Y1 <- get.Y(C1, C2, C3, C4, X = 1, eps)
  Y0 <- get.Y(C1, C2, C3, C4, X = 0, eps)
  Y <- get.Y(C1, C2, C3, C4, X, eps)

  if(do.pop){
    # return the pscore & counterfactual outcomes
    yay <- data.frame(cbind(pscore, Y1, Y0))
  } else {

    if(do.misp & my.transform){
      # BALZER & WESTLING transformation of the confounders
      C1 <- exp(C1 / 2)
      C2 <- C2 / (1 + exp(C1)) + 10
      C3 <- (C1 * C3 / 25 + 0.6) ^ 3
      C4 <- (C2 + C4 + 20) ^ 2

    } else if (do.misp & !my.transform){
      # NAIMI et al. transformation of the confounders
      Z1 <- C1
      Z2 <- C2
      Z3 <- C3
      Z4 <- C4

      C1 <- exp(Z1 / 2)
      C2 <- Z2 / (1 + exp(Z1)) + 10
      C3 <- (Z1 * Z3 / 25 + 0.6) ^ 3
      C4 <- (Z2 + Z4 + 20) ^ 2
    }
    # return observed data
    yay <- data.frame(cbind(C1, C2, C3, C4, X, Y))

  }
  yay
}
