
#' Estimate ATE using cross-fit procedure for TMLE estimator
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
#' sum(1:10)
#'
#'


tmle_single_p = function(data, exposure, outcome, covarsT, covarsO, learners, control, n_split){

    # Split sample

    splits_p <- sample(rep(1:n_split, diff(floor(nrow(data) * c(0:n_split/n_split)))))

    data_pp = data_p = data %>% mutate(s=splits_p)

    # Create nested dataset

    dat_nested_p <- data_p %>%
      group_by(s) %>%
      tidyr::nest()

    n_s = sapply(dat_nested_p$data, function(x) nrow(x))

    # P-score model

    pi_fitter <- function(df){
      SuperLearner::SuperLearner(Y=as.matrix(df[, exposure]), X=df[, covarsT], family=binomial(), SL.library=learners, cvControl=control)
    }

    dat_nested_p <- dat_nested_p %>%
      mutate(pi_fit=map(data, pi_fitter))


    # Calc p-scores using each split
    pi = H1 = H0 = list()
    k = dim(data_p)[2]
    for(i in 1:n_split){
      pi[[i]] = predict(dat_nested_p$pi_fit[[i]], newdata = data_pp[, covarsT])$pred
      pi[[i]] = ifelse(pi[[i]] < 0.025, 0.025, ifelse(pi[[i]] > 0.975, 0.975, pi[[i]]))
      H1[[i]] = data_pp[, exposure]/pi[[i]]
      H0[[i]] = (1 - data_pp[, exposure])/(1 - pi[[i]])
      data_p = suppressMessages(bind_cols(data_p, pi[[i]], H1[[i]], H0[[i]]))

    }

    names(data_p) <-   c(names(data_p)[1:k], paste0(rep(c("pi", "H1_", "H0_"), times = n_split), rep(1:n_split, each = 3)))

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
    mu_name = c(paste0(rep(c("mu_", "mu1_", "mu0_"), times = n_split), rep(1:n_split, each = 3)))
    mu = mu1 = mu0 = list()
    for(i in 1:n_split){
      mu[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = data_pp[, c(exposure, covarsO)])$pred
      mu[[i]] = ifelse(mu[[i]] == 0, 1e-17, ifelse(mu[[i]] == 1, 1-1e-17, mu[[i]]))
      mu1[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat1_p[, c(exposure, covarsO)])$pred
      mu1[[i]] = ifelse(mu1[[i]] == 0, 1e-17, ifelse(mu1[[i]] == 1, 1-1e-17, mu1[[i]]))
      mu0[[i]] = predict(dat_nested_p$mu_fit[[i]], newdata = dat0_p[, c(exposure, covarsO)])$pred
      mu0[[i]] = ifelse(mu0[[i]] == 0, 1e-17, ifelse(mu0[[i]] == 1, 1-1e-17, mu0[[i]]))
      data_p = suppressMessages(bind_cols(data_p, mu[[i]], mu1[[i]], mu0[[i]]))
    }

    names(data_p) = c(names(data_p)[1:k1], mu_name)


    pi_p <- data_p %>%
      select(s, paste0("pi", 1:n_split))

    mu_p <- data_p %>%
      select(s, paste0("mu_", 1:n_split))

    mu1_p <- data_p %>%
      select(s, paste0("mu1_", 1:n_split))

    mu0_p <- data_p %>%
      select(s, paste0("mu0_", 1:n_split))

    Y_p <-  data_p %>%
      select(s, Y)

    X_p <-  data_p %>%
      select(s, exposure)

    H_p <-  data_p %>%
      select(s, paste0(rep(c("H0_", "H1_"), times = n_split), rep(1:n_split, each = 2)))

    epsilon = mu0_1 = mu1_1 = iid = list()
    for(i in 1:n_split){
      iid[[i]] = sample(setdiff(1:n_split, i) , 2, replace = FALSE)

      h0 = pull(H_p %>% filter(s==i), paste0("H0_", iid[[i]][1]))
      h1 = pull(H_p %>% filter(s==i), paste0("H1_", iid[[i]][1]))
      muu = pull(mu_p %>% filter(s==i), paste0("mu_", iid[[i]][2]))

      epsilon[[i]] <- coef(glm(Y ~ -1 + h0 + h1 + offset(qlogis(muu)),
                               data = data_p %>% filter(s==i), family = binomial))

      mu0_1[[i]] = plogis(qlogis(pull(mu0_p, paste0("mu0_", iid[[i]][2]))) + epsilon[[i]][1] / (1 - pull(pi_p, paste0("pi", iid[[i]][1]))))
      mu1_1[[i]] = plogis(qlogis(pull(mu1_p, paste0("mu1_", iid[[i]][2]))) + epsilon[[i]][2] / pull(pi_p, paste0("pi", iid[[i]][1])))

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
        summarise(across(everything(), mean))

      r0[i] = data_p %>% filter(s==i) %>% select(paste0("mu0_1_", i)) %>%
        summarise(across(everything(), mean))

      rd[i] = r1[[i]] - r0[[i]]

      nm = pull(data_p, exposure)/pull(data_p, paste0("pi", iid[[i]][1]))
      dm = pull(data_p, Y) - pull(data_p, paste0("mu1_1_", iid[[i]][2]))
      ad = pull(data_p, paste0("mu1_1_", iid[[i]][2])) - r1[[i]]

      if1[[i]] = nm*dm + ad
      v1[i] = var(if1[[i]])/n_s[i]


      nm = (1 - pull(data_p, exposure))/(1 - pull(data_p, paste0("pi", iid[[i]][1])))
      dm = dplyr::pull(data_p, Y) - pull(data_p, paste0("mu0_1_", iid[[i]][2]))
      ad = dplyr::pull(data_p, paste0("mu0_1_", iid[[i]][2])) - r0[[i]]

      if0[[i]] = nm*dm + ad
      v0[i] = var(if0[[i]])/n_s[i]

      ifd[[i]] = if1[[i]] - if0[[i]]
      vd[i] = var(ifd[[i]])/n_s[i]

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
