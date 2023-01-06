#' Estimate Average Treatment Effect (ATE) using AIPW estimator with cross-fit algorithm (single repetition)
#'
#' @param data a data frame of tibble
#' @param exposure name of exposure variable
#' @param outcome name of outcome variable
#' @param covarsT a vector of names of covaraites for treatment model
#' @param covarsO a vector of names of covaraites for outcome model
#' @param family.y it is the family for outcome model. It can `binomial() (default)` or `"gaussian"`
#' @param learners similar as \code{SL.library()} in `SuperLearner` package.
#' @param control similar as  \code{cvControl()} in `SuperLearner` package.
#' @param n_split number of splits used, default `n_split = 3`
#' @param rand_split logical value; if be FALSE `(default)`, discordant splits for exposure and outcome model are chosen systematically; otherwise chosen randomly.
#' @param gbound value between (0,1) for truncation of predicted probabilities. See \code{tmle::tmle()} for more information.
#' @param alpha used to keep predicted initial values bounded away from (0,1) for logistic fluctuation.
#' @param seed numeric value to reproduce the splits distribution
#' @return a tibble of estimates.
#'
#' @import dplyr tibble tidyr purrr furrr parallel
#'
#' @importFrom stats binomial coef glm median plogis predict qlogis var
#'
#' @export
#'
#' @examples
#'
#' # See the README file for details
#'
#'
#' sum(1:5)
#'
aipw_single_p <- function(data,
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
                          seed){

  #suppressMessages(require(SuperLearner))
  # Split sample
  set.seed(seed)
  splits_p <- sample(rep(1:n_split, diff(floor(nrow(data) * c(0:n_split/n_split)))))

  data_p <- data %>% mutate(s=splits_p)

  # Create nested dataset

  dat_nested_p <- data_p %>%
    group_by(s) %>%
    nest()

  # P-score model

  pi_fitter <- function(df){
    SuperLearner::SuperLearner(Y=as.matrix(df[, exposure]),
                 X=df[, covarsT],
                 family=binomial(),
                 SL.library=learners,
                 cvControl=control)
  }

  dat_nested_p <- dat_nested_p %>%
    mutate(pi_fit=map(data, pi_fitter))


  # Calc p-scores using each split
  pi = list()
  k = dim(data_p)[2]
  for(i in 1:n_split){
    pi[[i]] = predict(dat_nested_p$pi_fit[[i]], newdata = data_p[, covarsT])$pred
    pi[[i]] = ifelse(pi[[i]] < gbound, gbound, ifelse(pi[[i]] > (1-gbound), (1-gbound), pi[[i]]))
    data_p = suppressMessages(bind_cols(data_p, pi[[i]]))

  }
  names(data_p) <- c(names(data_p)[1:k], paste0("pi", 1:n_split))

  #Outcome model
  if(family.y == "binomial"){
    mu_fitter <- function(df){
      SuperLearner::SuperLearner(Y=as.matrix(df[, outcome]),
                                 X=df[, c(exposure, covarsO)],
                                 family=binomial(),
                                 SL.library=learners,
                                 cvControl=control)
    }
  }

  if(family.y == "gaussian"){
    mu_fitter <- function(df){
      SuperLearner::SuperLearner(Y=as.matrix(df[, outcome]),
                                 X=df[, c(exposure, covarsO)],
                                 family="gaussian",
                                 SL.library=learners,
                                 cvControl=control)
    }
  }





  dat_nested_p <- dat_nested_p %>%
    mutate(mu_fit=map(data, mu_fitter))

  # Calc mu using each split
  dat1_p = dat0_p = data_p

  dat1_p[, exposure] = 1

  dat0_p[, exposure] = 0

  k1 = dim(data_p)[2]
  mu_name = c(paste0(rep(c("mu1_", "mu0_"), times = n_split), rep(1:n_split, each = 2)))
  mu1 = mu0 = list()

  if(family.y == "binomial"){

    for(i in 1:n_split){
      mu1[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat1_p %>% select(exposure, covarsO))$pred
      mu1[[i]] = ifelse(mu1[[i]] == 0, alpha, ifelse(mu1[[i]] == 1, 1-alpha, mu1[[i]]))
      mu0[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat0_p %>% select(exposure, covarsO))$pred
      mu0[[i]] = ifelse(mu0[[i]] == 0, alpha, ifelse(mu0[[i]] == 1, 1-alpha, mu0[[i]]))
      data_p = suppressMessages(bind_cols(data_p, mu1[[i]], mu0[[i]]))
    }
  }

  if(family.y == "gaussian"){

    for(i in 1:n_split){
      mu1[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat1_p %>% select(exposure, covarsO), type = "response")$pred
      mu1[[i]] = ifelse(mu1[[i]] == 0, alpha, ifelse(mu1[[i]] == 1, 1-alpha, mu1[[i]]))
      mu0[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat0_p %>% select(exposure, covarsO), type = "response")$pred
      mu0[[i]] = ifelse(mu0[[i]] == 0, alpha, ifelse(mu0[[i]] == 1, 1-alpha, mu0[[i]]))
      data_p = suppressMessages(bind_cols(data_p, mu1[[i]], mu0[[i]]))
    }
  }

  names(data_p) = c(names(data_p)[1:k1], mu_name)

  data_p <-  data_p %>% mutate(across(paste0("pi", 1:n_split), ~ pull(data_p, exposure)/.x +
                                        (1-pull(data_p, exposure))/(1-.x) , .names = "ipw_{.col}"))

  W_p <- data_p %>%
    select(s, paste0("ipw_pi", 1:n_split))

  mu1_p <- data_p %>%
    select(s, paste0("mu1_", 1:n_split))

  mu0_p <- data_p %>%
    select(s, paste0("mu0_", 1:n_split))

  Y_p <-  data_p %>%
    select(s, Y)

  X_p <-  data_p %>%
    select(s, exposure)

  a1 <- a0 <- list()
  if(rand_split == TRUE){

    for(i in 1:n_split){
      iid = sample(setdiff(1:n_split, i) , 2, replace = FALSE)
      a1[[i]] = unlist(X_p[X_p$s==i, -1] * W_p[W_p$s == i, paste0("ipw_pi", iid[1])] * (Y_p[Y_p$s==i, -1] - mu1_p[mu1_p$s == i, paste0("mu1_", iid[2])]) +
                         mu1_p[mu1_p$s == i, paste0("mu1_", iid[2])])
      a0[[i]] = unlist((1-X_p[X_p$s==i, -1]) * W_p[W_p$s == i, paste0("ipw_pi", iid[1])] * (Y_p[Y_p$s==i, -1] - mu0_p[mu0_p$s == i, paste0("mu0_", iid[2])]) +
                         mu0_p[mu0_p$s == i, paste0("mu0_", iid[2])])

    }

  }

  if(rand_split == FALSE){
    pi_id = c(2:n_split, 1); mu_id = c(3:n_split, 1, 2)
    for(i in 1:n_split){
      a1[[i]] = unlist(X_p[X_p$s==i, -1] * W_p[W_p$s == i, paste0("ipw_pi", pi_id[i])] * (Y_p[Y_p$s==i, -1] - mu1_p[mu1_p$s == i, paste0("mu1_", mu_id[i])]) +
                         mu1_p[mu1_p$s == i, paste0("mu1_", mu_id[i])])
      a0[[i]] = unlist((1-X_p[X_p$s==i, -1]) * W_p[W_p$s == i, paste0("ipw_pi", pi_id[i])] * (Y_p[Y_p$s==i, -1] - mu0_p[mu0_p$s == i, paste0("mu0_", mu_id[i])]) +
                         mu0_p[mu0_p$s == i, paste0("mu0_", mu_id[i])])

    }

  }




  r1 = sapply(a1, function(x) mean(x, na.rm = TRUE))
  r0 = sapply(a0, function(x) mean(x, na.rm = TRUE))
  rd = r1 - r0

  if_a1 = lapply(a1, function(x) x - mean(x, na.rm = TRUE))
  if_a0 = lapply(a0, function(x) x - mean(x, na.rm = TRUE))

  n_a1 = lapply(a1, function(x) length(x))
  rd1 = list()
  for(i in 1:n_split){
    rd1[[i]] = rep(rd[i], each = n_a1[[i]])
  }

  ifd = mapply(function(x, y, z) x - y - z, x = a1, y = a0, z = rd1, SIMPLIFY = FALSE)

  #Results
  r1_f = mean(r1)
  r0_f = mean(r0)
  rd_f = mean(rd)

  v1 = mean(sapply(if_a1, function(x) var(x))/unlist(n_a1))
  v0 = mean(sapply(if_a0, function(x) var(x))/unlist(n_a1))
  vd = mean(sapply(ifd, function(x) var(x))/unlist(n_a1))

  results <- tibble(r1=r1_f, r0=r0_f, rd=rd_f, v1=v1, v0=v0, vd=vd)
  return(results)
}


