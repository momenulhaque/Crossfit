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
                          family.y = "binomial",
                          learners,
                          control,
                          n_split = 5,
                          rand_split = FALSE,
                          gbound = 0.025,
                          alpha = 1e-17,
                          seed=146){

  # suppressMessages(require(SuperLearner))
  # Split sample

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

  #### step: 1.2.2 ############

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


  # Calculation of mu using each split

  dat1_p = dat0_p = data_pp

  dat1_p[, exposure] = 1

  dat0_p[, exposure] = 0

  k1 = dim(data_p)[2]

  mu_name = c(paste0(rep(c("mu_", "mu1_", "mu0_"), times = n_split), rep(1:n_split, each = 3)))
  mu = mu1 = mu0 = list()

  ########## 1.3.2 ##########################

  if(family.y == "binomial"){
    for(i in 1:n_split){
      mu[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = data_pp[, c(exposure, covarsO)], type = "response")$pred
      mu[[i]] = ifelse(mu[[i]] == 0, alpha, ifelse(mu[[i]] == 1, 1-alpha, mu[[i]]))
      mu1[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat1_p[, c(exposure, covarsO)], type = "response")$pred
      mu1[[i]] = ifelse(mu1[[i]] == 0, alpha, ifelse(mu1[[i]] == 1, 1-alpha, mu1[[i]]))
      mu0[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat0_p[, c(exposure, covarsO)], type = "response")$pred
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


  ################## Step: 1.3.3 #######################################

  epsilon = iid = list(); mu0_1 = mu1_1 = mu0_1s = mu1_1s = mu0_1ss = mu1_1ss = list()
  if(rand_split == TRUE){
    for(i in 1:n_split){
      iid[[i]] = sample(setdiff(1:n_split, i) , 2, replace = FALSE)

      h0 = H_p[ , paste0("H0_", iid[[i]][1])]

      h1 = H_p[ , paste0("H1_", iid[[i]][1])]

      muu = mu_p[ , paste0("mu_", iid[[i]][2])]

      y =  Y_p[ , outcome]

      dd <- data.frame(y=y, h0=h0, h1=h1, muu=muu)

      epsilon[[i]] <- coef(glm(dd[,1] ~ -1 + dd[,2] + dd[,3] + offset(qlogis(dd[,4])), family = binomial()))

      mu0_1s[[i]] = plogis(qlogis(pull(mu0_p, paste0("mu0_",  iid[[i]][2]))) + epsilon[[i]][1] / (1 - pull(pi_p, paste0("pi", iid[[i]][1]))))

      mu1_1s[[i]] = plogis(qlogis(pull(mu1_p, paste0("mu1_", iid[[i]][2]))) + epsilon[[i]][2] / pull(pi_p, paste0("pi", iid[[i]][1])))

    }

  }


  if(rand_split == FALSE){
    pi_id = c(2:n_split, 1); mu_id = c(3:n_split, 1, 2)
    for(i in 1:n_split){

      h0 = H_p[, paste0("H0_", pi_id[i])]

      h1 = H_p[, paste0("H1_", pi_id[i])]

      muu = mu_p[, paste0("mu_", mu_id[i])]

      y =  Y_p[, outcome]

      dd <- data.frame(y=y, h0=h0, h1=h1, muu=muu)

      epsilon[[i]] <- coef(glm(dd[,1] ~ -1 + dd[,2] + dd[,3] + offset(qlogis(dd[,4])), family = binomial()))

      mu0_1s[[i]] = plogis(qlogis(pull(mu0_p, paste0("mu0_",  mu_id[i]))) + epsilon[[i]][1] / (1 - pull(pi_p, paste0("pi", pi_id[i]))))

      mu1_1s[[i]] = plogis(qlogis(pull(mu1_p, paste0("mu1_", mu_id[i]))) + epsilon[[i]][2] / pull(pi_p, paste0("pi", pi_id[i])))

    }
  }

  mu0_1_p = suppressMessages(bind_cols(mu0_1s))
  mu1_1_p = suppressMessages(bind_cols(mu1_1s))

  dat_name = names(data_p)

  data_p = suppressMessages(bind_cols(data_p, mu0_1_p, mu1_1_p))

  names(data_p) = c(dat_name, paste0("mu0_1_", 1:n_split), paste0("mu1_1_", 1:n_split))



  r1 = r0 = rd = NULL
  r1_0 <- list()
  if1 = if0 = ifd = list()
  v1 = v0 = vd = NULL
  rdc <- NULL


  mu_summ <- list(); n <- NULL


  ############### step 1.3.4 - 1.3.5 ######################

  if(rand_split == FALSE){
    for(i in 1:n_split){
      data_p$ss <- data_p$s
      mu_summ[[i]] <- data_p %>% filter(s %in% i) %>%
        select(exposure, outcome, paste0("pi", i), paste0("mu1_1_", i), paste0("mu0_1_", i), ss)
      names(mu_summ[[i]]) <- c("x", "y", "pi", "mu1", "mu0", "s")
      n[i] <- nrow(mu_summ[[i]])
    }
  }


  if(rand_split == TRUE){

    for(i in 1:n_split){
      data_p$ss <- data_p$s
      mu_summ[[i]] <- data_p %>% filter(s %in% i) %>%
        select(exposure, outcome, paste0("pi", i), paste0("mu1_1_", i), paste0("mu0_1_", i), ss)
      names(mu_summ[[i]]) <- c("x", "y", "pi", "mu1", "mu0", "s")
      n[i] <- nrow(mu_summ[[i]])
    }
  }


  mu_sum <- bind_rows(mu_summ)
  mu_sum$ss <- rep(1:n_split, times = n)




  ############ step 1.6 ########################

  for(i in 1:n_split){
    rdc[i] <- mean(mu_sum[mu_sum$ss == i, ]$mu1) - mean(mu_sum[mu_sum$ss == i, ]$mu0)
    D1 = with(mu_sum[mu_sum$ss == i, ], (x/pi)*(y-mu1) + (mu1 - mean(mu1)))
    D0 = with(mu_sum[mu_sum$ss == i, ], ((1-x)/(1-pi))*(y-mu0) + (mu0 - mean(mu0)))
    EIC = D1-D0
    vd[i] <- var(EIC)/nrow(data)
  }


  rd = mean(rdc)
  var1 = mean(vd)

  res <- tibble(rd=rd, var = var1)


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

