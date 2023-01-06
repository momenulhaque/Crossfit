#' Estimate Average Treatment Effect (ATE) using TMLE estimator using cross-fit algorithm (single repetition)
#'
#' @param data a data frame of tibble
#' @param exposure name of exposure variable
#' @param outcome name of outcome variable
#' @param covarsT a vector of names of covaraites for treatment model
#' @param covarsO a vector of names of covaraites for outcome model
#' @param family.y it is the family for outcome model. It can `"binomial" (default)` or `"gaussian"`
#' @param learners similar as \code{SL.library()} in `SuperLearner` package.
#' @param control similar as  \code{cvControl()} in `SuperLearner` package.
#' @param n_split number of splits used, default `n_split = 3`
#' @param rand_split logical value; if be FALSE `(default)`, discordant splits for exposure and outcome model are chosen systematically; otherwise chosen randomly.
#' @param gbound value between (0,1) for truncation of predicted probabilities. See \code{tmle::tmle()} for more information.
#' @param alpha used to keep predicted initial values bounded away from (0,1) for logistic fluctuation.
#' @param seed numeric value to reproduce the splits distribution
#' @return a tibble of estimates.
#'
#' @import dplyr tibble tidyr purrr furrr tmle
#'
#' @importFrom stats binomial coef glm median plogis predict qlogis var
#'
#'
#'
#' @export
#'
#' @examples
#'
#' # See the README file for details
#'
#' sum(1:5)
#'
#'
#'
#'
tmle_single_p = function(data,
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

 # suppressMessages(require(SuperLearner))
  # Split sample
  set.seed(seed)
  splits_p <- sample(rep(1:n_split, diff(floor(nrow(data) * c(0:n_split/n_split)))))

  data_pp = data_p = data %>% mutate(s=splits_p)

  # Create nested dataset

  dat_nested_p <- data_p %>%
    group_by(s) %>%
    tidyr::nest()

  n_s = sapply(dat_nested_p$data, function(x) nrow(x))

  # P-score model

  pi_fitter <- function(df){
    SuperLearner::SuperLearner(Y=as.matrix(df[, exposure]),
                               X=df[, covarsT],
                               family=binomial(),
                               SL.library=learners,
                               cvControl=control)
  }

  dat_nested_p <- dat_nested_p %>% mutate(pi_fit = map(data, pi_fitter))


  # Calc p-scores using each split
  pi = H1 = H0 = list()
  k = dim(data_p)[2]
  for(i in 1:n_split){
    pi[[i]] = predict(dat_nested_p$pi_fit[[i]], newdata = data_pp[, covarsT])$pred
    pi[[i]] = ifelse(pi[[i]] < gbound, gbound, ifelse(pi[[i]] > (1-gbound), (1-gbound), pi[[i]]))
    H1[[i]] = data_pp[, exposure]/pi[[i]]
    H0[[i]] = (1 - data_pp[, exposure])/(1 - pi[[i]])
    data_p = suppressMessages(bind_cols(data_p, pi[[i]], H1[[i]], H0[[i]]))

  }

  names(data_p) <-   c(names(data_p)[1:k], paste0(rep(c("pi", "H1_", "H0_"),
                                                      times = n_split), rep(1:n_split, each = 3)))

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
  mu_name = c(paste0(rep(c("mu_", "mu1_", "mu0_"), times = n_split), rep(1:n_split, each = 3)))
  mu = mu1 = mu0 = list()

  if(family.y == "binomial"){
      for(i in 1:n_split){
        mu[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = data_pp[, c(exposure, covarsO)])$pred
        mu[[i]] = ifelse(mu[[i]] == 0, alpha, ifelse(mu[[i]] == 1, 1-alpha, mu[[i]]))
        mu1[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat1_p[, c(exposure, covarsO)])$pred
        mu1[[i]] = ifelse(mu1[[i]] == 0, alpha, ifelse(mu1[[i]] == 1, 1-alpha, mu1[[i]]))
        mu0[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat0_p[, c(exposure, covarsO)])$pred
        mu0[[i]] = ifelse(mu0[[i]] == 0, alpha, ifelse(mu0[[i]] == 1, 1-alpha, mu0[[i]]))
        data_p = suppressMessages(bind_cols(data_p, mu[[i]], mu1[[i]], mu0[[i]]))
      }
  }
  if(family.y == "gaussian"){
    for(i in 1:n_split){
      mu[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = data_pp[, c(exposure, covarsO)], type = "response")$pred
      mu[[i]] = ifelse(mu[[i]] <= 0, alpha, ifelse(mu[[i]] >= 1, 1-alpha, mu[[i]]))
      mu1[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat1_p[, c(exposure, covarsO)], type = "response")$pred
      mu1[[i]] = ifelse(mu1[[i]] <= 0, alpha, ifelse(mu1[[i]] >= 1, 1-alpha, mu1[[i]]))
      mu0[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat0_p[, c(exposure, covarsO)], type = "response")$pred
      mu0[[i]] = ifelse(mu0[[i]] <= 0, alpha, ifelse(mu0[[i]] >= 1, 1-alpha, mu0[[i]]))
      data_p = suppressMessages(bind_cols(data_p, mu[[i]], mu1[[i]], mu0[[i]]))
    }
  }


  names(data_p) = c(names(data_p)[1:k1], mu_name)


  pi_p <- suppressWarnings(data_p %>%
                             select(s, paste0("pi", 1:n_split)))

  mu_p <- suppressWarnings(data_p %>%
                             select(s, paste0("mu_", 1:n_split)))

  mu1_p <- suppressWarnings(data_p %>%
                              select(s, paste0("mu1_", 1:n_split)))

  mu0_p <- suppressWarnings(data_p %>%
                              select(s, paste0("mu0_", 1:n_split)))

  Y_p <-  suppressWarnings(data_p %>%
                             select(s, outcome))


  X_p <- suppressWarnings(data_p %>%
                            select(s, exposure))

  H_p <-  data_p %>%
    select(s, paste0(rep(c("H0_", "H1_"), times = n_split), rep(1:n_split, each = 2)))

  epsilon = mu0_1 = mu1_1 = iid = list()


  if(rand_split == TRUE){
    for(i in 1:n_split){
      iid[[i]] = sample(setdiff(1:n_split, i) , 2, replace = FALSE)

      h0 = pull(H_p %>% filter(s==i), paste0("H0_", iid[[i]][1]))
      h1 = pull(H_p %>% filter(s==i), paste0("H1_", iid[[i]][1]))
      muu = pull(mu_p %>% filter(s==i), paste0("mu_", iid[[i]][2]))
      y = pull(data_p %>% filter(s==i), outcome)

#      if(family.y == "binomial"){
        epsilon[[i]] <- coef(glm(y ~ -1 + h0 + h1 + offset(qlogis(muu)),
                                 data = data_p %>% filter(s==i), family = binomial()))
        mu0_1[[i]] = plogis(qlogis(pull(mu0_p, paste0("mu0_", iid[[i]][2]))) + epsilon[[i]][1] / (1 - pull(pi_p, paste0("pi", iid[[i]][1]))))
        mu1_1[[i]] = plogis(qlogis(pull(mu1_p, paste0("mu1_", iid[[i]][2]))) + epsilon[[i]][2] / pull(pi_p, paste0("pi", iid[[i]][1])))

    }
  }

  if(rand_split == FALSE){
    pi_id = c(2:n_split, 1); mu_id = c(3:n_split, 1, 2)
    for(i in 1:n_split){

      h0 = pull(H_p %>% filter(s==i), paste0("H0_", pi_id[i]))
      h1 = pull(H_p %>% filter(s==i), paste0("H1_", pi_id[i]))
      muu = pull(mu_p %>% filter(s==i), paste0("mu_", mu_id[i]))
      y = pull(data_p %>% filter(s==i), outcome)

        epsilon[[i]] <- coef(glm(y ~ -1 + h0 + h1 + offset(qlogis(muu)),
                                 data = data_p %>% filter(s==i), family = binomial()))
        mu0_1[[i]] = plogis(qlogis(pull(mu0_p, paste0("mu0_", mu_id[i]))) + epsilon[[i]][1] / (1 - pull(pi_p, paste0("pi", pi_id[i]))))
        mu1_1[[i]] = plogis(qlogis(pull(mu1_p, paste0("mu1_", mu_id[i]))) + epsilon[[i]][2] / pull(pi_p, paste0("pi", pi_id[i])))


    }
  }





  mu0_1_p = suppressMessages(bind_cols(mu0_1))
  mu1_1_p = suppressMessages(bind_cols(mu1_1))

  dat_name = names(data_p)

  data_p = suppressMessages(bind_cols(data_p, mu0_1_p, mu1_1_p))

  names(data_p) = c(dat_name, paste0("mu0_1_", 1:n_split), paste0("mu1_1_", 1:n_split))

  r1 = r0 = rd = NULL
  if1 = if0 = ifd = list()
  v1 = v0 = vd = NULL
  for(i in 1:n_split){
    r1[i] = data_p %>% filter(s==i) %>% select(paste0("mu1_1_", i)) %>%
      summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

    r0[i] = data_p %>% filter(s==i) %>% select(paste0("mu0_1_", i)) %>%
      summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

    rd[i] = r1[[i]] - r0[[i]]

    if(rand_split == TRUE){
      nm = pull(data_p, exposure)/pull(data_p, paste0("pi", iid[[i]][1]))
      dm = pull(data_p, outcome) - pull(data_p, paste0("mu1_1_", iid[[i]][2]))
      ad = pull(data_p, paste0("mu1_1_", iid[[i]][2])) - r1[[i]]
    }


    if(rand_split == FALSE){
      nm = pull(data_p, exposure)/pull(data_p, paste0("pi", pi_id[i]))
      dm = pull(data_p, outcome) - pull(data_p, paste0("mu1_1_", mu_id[i]))
      ad = pull(data_p, paste0("mu1_1_", mu_id[i])) - r1[[i]]
    }



    if1[[i]] = nm*dm + ad

    v1[i] = var(if1[[i]], na.rm = TRUE)/n_s[i]

    if(rand_split == TRUE){
      nm = (1 - pull(data_p, exposure))/(1 - pull(data_p, paste0("pi", iid[[i]][1])))
      dm = dplyr::pull(data_p, outcome) - pull(data_p, paste0("mu0_1_", iid[[i]][2]))
      ad = dplyr::pull(data_p, paste0("mu0_1_", iid[[i]][2])) - r0[[i]]
    }

    if(rand_split == FALSE){
      nm = (1 - pull(data_p, exposure))/(1 - pull(data_p, paste0("pi", pi_id[i])))
      dm = dplyr::pull(data_p, outcome) - pull(data_p, paste0("mu0_1_", mu_id[i]))
      ad = dplyr::pull(data_p, paste0("mu0_1_", mu_id[i])) - r0[[i]]
    }


    if0[[i]] = nm*dm + ad
    v0[i] = var(if0[[i]], na.rm = TRUE)/n_s[i]

    ifd[[i]] = if1[[i]] - if0[[i]]
    vd[i] = var(ifd[[i]], na.rm = TRUE)/n_s[i]

  }

  r1 = mean(unlist(r1))
  r0 = mean(unlist(r0))
  rd = mean(unlist(rd))

  v1 = mean(v1)
  v0 = mean(v0)
  vd = mean(vd)

  results <- tibble(r1, r0, rd, v1, v0, vd)
  return(results)
}

environment(tmle_single_p) <- .GlobalEnv
