


#' Estimate ATE using cross-fit procedure for AIPW estimator
#'
#' @param data a data frame of tibble
#' @param exposure name of exposure variable
#' @param outcome name of outcome variable
#' @param covarsT a vector of names of covaraites for treatment model
#' @param covarsO a vector of names of covaraites for outcome model
#' @param learners similar as \code{\link[Superlearner:SL.library()]{Superlearner::SL.library()}}
#' @param control similar as  \code{\link[Superlearner:cvControl()]{Superlearner::cvControl()}}
#' @param n_split number of splits
#'
#' @return a tibble of estimates
#'
#' @export
#'
#' @examples
#'
#' sum(1:5)
#'
aipw_single_p <- function(data, exposure, outcome, covarsT, covarsO, learners, control, n_split){

  # Split sample

  splits_p <- sample(rep(1:n_split, diff(floor(nrow(data) * c(0:n_split/n_split)))))

  data_p <- data %>% mutate(s=splits_p)

  # Create nested dataset

  dat_nested_p <- data_p %>%
    group_by(s) %>%
    nest()

  # P-score model

  pi_fitter <- function(df){
    SuperLearner::SuperLearner(Y=as.matrix(df[, exposure]), X=df[, covarsT], family=binomial(), SL.library=learners, cvControl=control)
  }

  dat_nested_p <- dat_nested_p %>%
    mutate(pi_fit=map(data, pi_fitter))


  # Calc p-scores using each split
  pi = list()
  k = dim(data_p)[2]
  for(i in 1:n_split){
    pi[[i]] = predict(dat_nested_p$pi_fit[[i]], newdata = data_p[, covarsT])$pred
    pi[[i]] = ifelse(pi[[i]] < 0.025, 0.025, ifelse(pi[[i]] > 0.975, 0.975, pi[[i]]))
    data_p = suppressMessages(bind_cols(data_p, pi[[i]]))

  }
  names(data_p) <- c(names(data_p)[1:k], paste0("pi", 1:n_split))

  #Outcome model
  mu_fitter <- function(df){
    SuperLearner::SuperLearner(Y=as.matrix(df[, outcome]), X=df[, c(exposure, covarsO)], family=binomial(), SL.library=learners, cvControl=control)
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
  for(i in 1:n_split){
    mu1[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat1_p %>% select(exposure, covarsO))$pred
    mu0[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat0_p %>% select(exposure, covarsO))$pred
    data_p = suppressMessages(bind_cols(data_p, mu1[[i]], mu0[[i]]))
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
  for(i in 1:n_split){
    iid = sample(setdiff(1:n_split, i) , 2, replace = FALSE)
    a1[[i]] = unlist(X_p[X_p$s==i, -1] * W_p[W_p$s == i, paste0("ipw_pi", iid[1])] * (Y_p[Y_p$s==i, -1] - mu1_p[mu1_p$s == i, paste0("mu1_", iid[2])]) +
                       mu1_p[mu1_p$s == i, paste0("mu1_", iid[2])])
    a0[[i]] = unlist((1-X_p[X_p$s==i, -1]) * W_p[W_p$s == i, paste0("ipw_pi", iid[1])] * (Y_p[Y_p$s==i, -1] - mu0_p[mu0_p$s == i, paste0("mu0_", iid[2])]) +
                       mu0_p[mu0_p$s == i, paste0("mu0_", iid[2])])

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

