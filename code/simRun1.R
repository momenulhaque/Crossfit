require(parallel)
require(plyr)

#load(file = "data/wrapperCode.RData")
numCores <- detectCores()
numCoresUsed <- numCores

cl <- makeCluster(numCoresUsed)

clusterEvalQ(cl,{
  #source the simulation code
  source("code/estimator.R", local = TRUE)
  source("code/estimator_p.R", local = TRUE)
  load(file = "data/wrapperCode.RData")
  
})

statin_sim_data <- read.csv("data/statin_sim_data.csv")

sims <- statin_sim_data %>%
  group_by(sim_id) %>%
  nest()

n = 2

aipw_result = clusterMap(cl, function(df, seed) aipw_sim(df, num_cf=3, seed), df=sims$data[1:n], seed=sims$sim_id[1:n])
saveRDS(aipw_result, file = "data/aipw_result.rds")

tmle_result = clusterMap(cl, function(df, seed) tmle_sim(df, num_cf=10, seed), df=sims$data[1:n], seed=sims$sim_id[1:n])
saveRDS(tmle_result, file = "data/tmle_result.rds")

aipw_result_p = clusterMap(cl, function(df, seed) aipw_sim_p(df, n_split = 3, num_cf=10, seed), df=sims$data[1:n], seed=sims$sim_id[1:n])
saveRDS(aipw_result, file = "data/aipw_result_p.rds")

tmle_result_p = clusterMap(cl, function(df, seed) tmle_sim_p(df, n_split = 3, num_cf=10, seed), df=sims$data[1:n], seed=sims$sim_id[1:n])
saveRDS(tmle_result, file = "data/tmle_result_p.rds")


aipw_result_p10 = clusterMap(cl, function(df, seed) aipw_sim_p(df, n_split = 10, num_cf=10, seed), df=sims$data[1:n], seed=sims$sim_id[1:n])
saveRDS(aipw_result_p10, file = "data/aipw_result_p10.rds")

tmle_result_p10 = clusterMap(cl, function(df, seed) tmle_sim_p(df, n_split = 10, num_cf=10, seed), df=sims$data[1:n], seed=sims$sim_id[1:n])
saveRDS(aipw_result_p10, file = "data/tmle_result_p10.rds")

stopCluster(cl)

aipw_resultF <- apply(dplyr::bind_rows(aipw_result), 2, mean)
tmle_resultF <- apply(dplyr::bind_rows(tmle_result), 2, mean)
aipw_resultF_p <- apply(dplyr::bind_rows(aipw_result_p), 2, mean)
tmle_resultF_p <- apply(dplyr::bind_rows(tmle_result_p), 2, mean)
aipw_resultF_p10 <- apply(dplyr::bind_rows(aipw_result_p10), 2, mean)
tmle_resultF_p10 <- apply(dplyr::bind_rows(tmle_result_p10), 2, mean)


result = bind_rows(aipw_resultF, aipw_resultF_p, aipw_resultF_p10, tmle_resultF, tmle_resultF_p,tmle_resultF_p10)
result$method <- c("aipw", "aipw_p", "aipw_p10", "tmle", "tmle_p", "tmle_p10")
saveRDS(result, file = "data/simResultF.Rds")


