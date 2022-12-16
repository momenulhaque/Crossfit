library(tidyverse)
library(tmle)
require(purrr)
require(furrr)
require(dplyr)
require(tidyr)
require(tibble)
require(SuperLearner)

aipw_single_p <- function(data, exposure, outcome, covarsT, covarsO, learners, control, n_split, rand_split = TRUE, seed = 145){
  
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
    SuperLearner(Y=as.matrix(df[, exposure]), X=df[, covarsT], family=binomial(), SL.library=learners, cvControl=control)
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
    select(s, outcome)
  
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



aipw_multiple_p <-function(data, exposure, outcome, covarsT, covarsO, learners, control, num_cf, n_split, rand_split = TRUE, seed = 145){
  
  #Initialize results
  runs <- tibble(r1=double(), r0=double(), rd=double(), v1=double(), v0=double(), vd=double())
  
  #Run on num_cf splits
  set.seed(seed)
  cf_seed = sample(num_cf)
  for(cf in 1:num_cf){
    seed = cf_seed[cf]
    runs <- bind_rows(runs, aipw_single_p(data, exposure, outcome, covarsT, covarsO, learners, control, n_split, rand_split, seed))
    
  }
  #Medians of splits
  medians <- apply(runs, 2, median)
  
  
  #Corrected variance terms
  runs <- runs %>%
    mutate(mv1 = v1 + (r1-medians[1])^2,
           mv0 = v0 + (r0-medians[2])^2,
           mvd = vd + (rd-medians[3])^2)
  
  results <- apply(runs, 2, median)
  return(results)
  
}




data = df


tmle_single_p = function(data, exposure, outcome, covarsT, covarsO, learners, control, n_split, rand_split = TRUE, seed = 145){
  
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
    select(s, outcome)
  
  X_p <-  data_p %>%
    select(s, exposure)
  
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
      
      epsilon[[i]] <- coef(glm(y ~ -1 + h0 + h1 + offset(qlogis(muu)),
                               data = data_p %>% filter(s==i), family = binomial))
      
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
      
      epsilon[[i]] <- coef(glm(outcome ~ -1 + h0 + h1 + offset(qlogis(muu)),
                               data = data_p %>% filter(s==i), family = binomial))
      
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
      summarise(across(everything(), mean))
    
    r0[i] = data_p %>% filter(s==i) %>% select(paste0("mu0_1_", i)) %>%
      summarise(across(everything(), mean))
    
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
    v1[i] = var(if1[[i]])/n_s[i]
    
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



tmle_multiple_p <- function(data, exposure, outcome, covarsT, covarsO, learners, control, num_cf, n_split, rand_split = TRUE, seed = 145){
  
  #Initialize results
  runs <- tibble(r1=double(), r0=double(), rd=double(), v1=double(), v0=double(), vd=double())
  
  #Run on num_cf splits  
  set.seed(seed)
  cf_seed = sample(num_cf)
  for(cf in 1:num_cf){
    seed = cf_seed[cf]
    runs <- bind_rows(runs, tmle_single_p(data, exposure, outcome, covarsT, covarsO, learners, control, n_split, rand_split, seed))
    
  }
  #Medians of splits
  medians <- apply(runs, 2, median)
  
  
  #Corrected variance terms
  runs <- runs %>%
    mutate(mv1 = v1 + (r1-medians[1])^2,
           mv0 = v0 + (r0-medians[2])^2,
           mvd = vd + (rd-medians[3])^2)
  
  results <- apply(runs, 2, median)
  return(results)
  
}



