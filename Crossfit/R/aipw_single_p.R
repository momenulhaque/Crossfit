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
#' @import dplyr tibble tidyr purrr furrr
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
                            family.y = "binomial",
                            learners,
                            control,
                            n_split = 5,
                            rand_split = FALSE,
                            gbound = 0.025,
                            alpha = 1e-17,
                            seed=146){

  #### step: 1.1 ############

  set.seed(seed)
  splits_p <- sample(rep(1:n_split, diff(floor(nrow(data) * c(0:n_split/n_split)))))

  data_pp = data_p = data %>% mutate(s=splits_p)  %>% arrange(s)

  # Create nested dataset

  dat_nested_p <- data_p %>%
    group_by(s) %>%
    tidyr::nest()


  #### step: 1.2.1 ############

  # P-score model

  pi_fitter <- function(df){
    SuperLearner::SuperLearner(Y=as.matrix(df[, exposure]),
                               X=df[, covarsT],
                               family=binomial(),
                               SL.library=learners,
                               cvControl=control)
  }

  dat_nested_p <- dat_nested_p %>% mutate(pi_fit = map(data, pi_fitter))



  ########## 1.3.1 ##########################
  pi = list()
  k = dim(data_p)[2]
  for(i in 1:n_split){
    pi[[i]] = predict(dat_nested_p$pi_fit[[i]], newdata = data_p[, covarsT])$pred
    pi[[i]] = ifelse(pi[[i]] < gbound, gbound, ifelse(pi[[i]] > (1-gbound), (1-gbound), pi[[i]]))
    data_p = suppressMessages(bind_cols(data_p, pi[[i]]))

  }
  names(data_p) <- c(names(data_p)[1:k], paste0("pi", 1:n_split))

  # Outcome model

  ############### step: 1.2.2 #############################

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





  dat_nested_p <- dat_nested_p %>% mutate(mu_fit=map(data, mu_fitter))

  # Calc mu using each split
  dat1_p = dat0_p = data_p

  dat1_p[, exposure] = 1

  dat0_p[, exposure] = 0

  k1 = dim(data_p)[2]
  mu_name = c(paste0(rep(c("mu1_", "mu0_"), times = n_split), rep(1:n_split, each = 2)))
  mu1 = mu0 = list()

  ########## 1.3.2 ##########################

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



  a1s <- a0s <- iid <- list()

  ################## Step: 1.3.3 #######################################

  if(rand_split == TRUE){

      for(i in 1:n_split){

        iid[[i]] = sample(setdiff(1:n_split, i), 2, replace = FALSE)

        a1s[[i]] = unlist(X_p[, exposure] * W_p[, paste0("ipw_pi", iid[[i]][1])] * (Y_p[, outcome] - mu1_p[, paste0("mu1_", iid[[i]][2])]) +
                           mu1_p[, paste0("mu1_", iid[[i]][2])])
        a0s[[i]] = unlist((1-X_p[, exposure]) * W_p[, paste0("ipw_pi", iid[[i]][1])] * (Y_p[, outcome] - mu0_p[, paste0("mu0_", iid[[i]][2])]) +
                           mu0_p[, paste0("mu0_", iid[[i]][2])])

      }

    }

    if(rand_split == FALSE){
      pi_id = c(2:n_split, 1); mu_id = c(3:n_split, 1, 2)
      for(i in 1:n_split){
        a1s[[i]] = unlist(X_p[, exposure] * W_p[, paste0("ipw_pi", pi_id[i])] * (Y_p[, outcome] - mu1_p[, paste0("mu1_", mu_id[i])]) +
                           mu1_p[, paste0("mu1_", mu_id[i])])
        a0s[[i]] = unlist((1-X_p[, exposure]) * W_p[, paste0("ipw_pi", pi_id[i])] * (Y_p[, outcome] - mu0_p[, paste0("mu0_", mu_id[i])]) +
                           mu0_p[, paste0("mu0_", mu_id[i])])

      }

    }


    a1 <- suppressMessages(bind_cols(a1s))
    a0 <- suppressMessages(bind_cols(a0s))

    dat_name = names(data_p)

    data_p = suppressMessages(bind_cols(data_p, a1, a0))
    names(data_p) = c(dat_name, paste0("a1_", 1:n_split), paste0("a0_", 1:n_split))

    rd = var.rd = NULL

    ############ step 1.6 and 1.7 ########################


    for(i in 1:n_split){

        a11 = data_p %>% filter(s==i) %>% select(paste0("a1_", i)) %>% pull()
        a10 = data_p %>% filter(s==i) %>% select(paste0("a0_", i)) %>% pull()

        rd[i] = mean(a11) - mean(a10)

        ifd = a11 - a10 - rd[i]

        var.rd[i] = var(ifd)/nrow(data)


    }


      rd = mean(rd)

      vd = mean(var.rd)

      res <- tibble(rd=rd, var=vd)


      df = dat_nested_p
      x_weight = bind_rows(lapply(df$pi_fit, function(x) x$coef))
      y_weight = bind_rows(lapply(df$mu_fit, function(x) x$coef))
      weight = bind_rows(x_weight, y_weight)
      weight$model = rep(c("x", "y"), each = n_split)
      weight$split = rep(1:n_split, times = 2)

      x_prev = unlist(bind_rows(lapply(df$data, function(x) colMeans(x[, exposure], na.rm = TRUE))))
      y_prev = unlist(bind_rows(lapply(df$data, function(x) colMeans(x[, outcome], na.rm = TRUE))))
      prev = c(x_prev, y_prev)
      weight$prev = prev

      fit = list(results=res, weight=weight)
      fit

  }
