#' Estimate Average Treatment Effect (ATE) using TMLE estimator using cross-fit algorithm for single repetition (generalization 2)
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
#' @param gbounds value between (0,1) for truncation of predicted probabilities of exposure. The defaults are 5/sqrt(n)log(n) and 1-5/sqrt(n)log(n). See \code{tmle::tmle()} for more information.
#' @param Qbounds used to keep predicted probabilities of outcomes values bounded away from (0,1). The defaults are 5e-04 and 1-5e-04.
#' @param seed numeric value to reproduce the splits distribution
#' @return A tibble containing risk (or mean) difference (`rd`), variance of the estimated `rd`.
#'
#' @import dplyr tibble tidyr purrr furrr tmle
#'
#' @importFrom stats binomial coef glm median plogis predict qlogis var
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
tmle_single_g2_p = function(data,
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
                            seed=146){

   #### step: 1.1 ############

  set.seed(seed)
  splits_p <- sample(rep(1:n_split, diff(floor(nrow(data) * c(0:n_split/n_split)))))

  if(family.y == "gaussian"){
    original.Y = outcome
    min.Y = min(data[, outcome])
    max.Y = max(data[, outcome])
    data$Y.star = (data[, outcome] - min.Y)/(max.Y-min.Y)
    outcome = "Y.star"
  }

  data_pp = data_p = data %>% mutate(s=splits_p)  %>% arrange(s)

  # Create nested dataset

  dat_nested_p <- data_p %>%
    group_by(s) %>%
    tidyr::nest()

  #### step: 1.2.1 ############

  # P-score model

  ##########################################
  pi_fitter <- function(df){
    SuperLearner::SuperLearner(Y=as.matrix(df[, exposure]),
                               X=df[, covarsT],
                               family=binomial(),
                               SL.library=learners,
                               cvControl=control)
  }


  ###############################################

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




  #################################################
  #### step: 1.2.1 ############

  dat_pi <- dat_mu <- pi_fitt <- mu_fitt <- list()

  if(rand_split == FALSE){

    rd_id <- NULL
    pi_id = mu_id = list()
    id = 1:n_split
    dif = (n_split-1)/2

    for(i in 1:n_split){
      rd_id[i] = i
      id_temp = c(id[id > rd_id[i]], id[id < rd_id[i]])
      pi_id[[i]] = id_temp[1:dif]
      mu_id[[i]] = setdiff(id_temp, pi_id[[i]])
    }


    for(i in 1:n_split){

      dat_pi[[i]] = data_p %>% filter(s %in% pi_id[[i]])
      pi_fitt[[i]] = pi_fitter(dat_pi[[i]])

      dat_mu[[i]] = data_p %>% filter(s %in% mu_id[[i]])
      mu_fitt[[i]] = mu_fitter(dat_mu[[i]])


    }
    dat_nested_p$pi_fit = pi_fitt
    dat_nested_p$mu_fit = mu_fitt
  }


  if(rand_split == TRUE){
    rd_id <- NULL;  pi_id = mu_id = list()
    for(i in 1:n_split){
      rd_id[i] <- i
      pi_id[[i]] = sample(setdiff(1:n_split, rd_id[i]) , (n_split-1)/2, replace = FALSE)
      mu_id[[i]] = setdiff(1:n_split, c(rd_id[i], pi_id[[i]]))
    }

    for(i in 1:n_split){

      dat_pi[[i]] = data_p %>% filter(s %in% pi_id[[i]])
      pi_fitt[[i]] = pi_fitter(dat_pi[[i]])

      dat_mu[[i]] = data_p %>% filter(s %in% mu_id[[i]])
      mu_fitt[[i]] = mu_fitter(dat_mu[[i]])


    }
    dat_nested_p$pi_fit = pi_fitt
    dat_nested_p$mu_fit = mu_fitt

  }







  ########## 1.3.1 ##########################
  pi = H1 = H0 = list()
  k = dim(data_p)[2]
  for(i in 1:n_split){
    pi[[i]] = predict(dat_nested_p$pi_fit[[i]], newdata = data_pp[, covarsT])$pred
    pi[[i]] = ifelse(pi[[i]] < gbounds, gbounds, ifelse(pi[[i]] > (1-gbounds), (1-gbounds), pi[[i]]))
    H1[[i]] = data_pp[, exposure]/pi[[i]]
    H0[[i]] = (1 - data_pp[, exposure])/(1 - pi[[i]])
    data_p = suppressMessages(bind_cols(data_p, pi[[i]], H1[[i]], H0[[i]]))

  }

  names(data_p) <-   c(names(data_p)[1:k], paste0(rep(c("pi", "H1_", "H0_"),
                                                      times = n_split), rep(1:n_split, each = 3)))


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
      mu[[i]] = ifelse(mu[[i]] < Qbounds, Qbounds, ifelse(mu[[i]] > 1-Qbounds, 1-Qbounds, mu[[i]]))
      mu1[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat1_p[, c(exposure, covarsO)], type = "response")$pred
      mu1[[i]] = ifelse(mu1[[i]] < Qbounds, Qbounds, ifelse(mu1[[i]] > 1-Qbounds, 1-Qbounds, mu1[[i]]))
      mu0[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat0_p[, c(exposure, covarsO)], type = "response")$pred
      mu0[[i]] = ifelse(mu0[[i]] < Qbounds, Qbounds, ifelse(mu0[[i]] > 1-Qbounds, 1-Qbounds, mu0[[i]]))
      data_p = suppressMessages(bind_cols(data_p, mu[[i]], mu1[[i]], mu0[[i]]))
    }
  }
  if(family.y == "gaussian"){
    for(i in 1:n_split){
      mu[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = data_pp[, c(exposure, covarsO)], type = "response")$pred
      mu[[i]] = ifelse(mu[[i]] < Qbounds, Qbounds, ifelse(mu[[i]] > 1-Qbounds, 1-Qbounds, mu[[i]]))
      mu1[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat1_p[, c(exposure, covarsO)], type = "response")$pred
      mu1[[i]] = ifelse(mu1[[i]] < Qbounds, Qbounds, ifelse(mu1[[i]] > 1-Qbounds, 1-Qbounds, mu1[[i]]))
      mu0[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat0_p[, c(exposure, covarsO)], type = "response")$pred
      mu0[[i]] = ifelse(mu0[[i]] < Qbounds, Qbounds, ifelse(mu0[[i]] > 1-Qbounds, 1-Qbounds, mu0[[i]]))
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

  epsilon = list(); mu0_1s = mu1_1s = list()

  for(i in 1:n_split){

    h0 = H_p[, paste0("H0_", i)]

    h1 = H_p[, paste0("H1_", i)]

    muu = mu_p[, paste0("mu_", i)]

    y =  Y_p[, outcome]

    dd <- data.frame(y=y, h0=h0, h1=h1, muu=muu)

    epsilon[[i]] <- coef(glm(dd[,1] ~ -1 + dd[,2] + dd[,3] + offset(qlogis(dd[,4])), family = binomial()))

    mu0_1s[[i]] = plogis(qlogis(pull(mu0_p, paste0("mu0_",  i))) + epsilon[[i]][1] / (1 - pull(pi_p, paste0("pi", i))))

    mu1_1s[[i]] = plogis(qlogis(pull(mu1_p, paste0("mu1_", i))) + epsilon[[i]][2] / pull(pi_p, paste0("pi", i)))

  }

  mu0_1_p = suppressMessages(bind_cols(mu0_1s))
  mu1_1_p = suppressMessages(bind_cols(mu1_1s))

  dat_name = names(data_p)

  data_p = suppressMessages(bind_cols(data_p, mu0_1_p, mu1_1_p))

  names(data_p) = c(dat_name, paste0("mu0_1_", 1:n_split), paste0("mu1_1_", 1:n_split))


  if(family.y == "gaussian"){
    for(i in 1:n_split){
      data_p[ ,paste0("mu1_1_", i)] <- data_p[ ,paste0("mu1_1_", i)]*(max.Y-min.Y) + min.Y
      data_p[ ,paste0("mu0_1_", i)] <- data_p[ ,paste0("mu0_1_", i)]*(max.Y-min.Y) + min.Y
    }
    outcome = original.Y
  }



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

  res <- data.frame(rd=rd, var = var1)

  return(res)

}

