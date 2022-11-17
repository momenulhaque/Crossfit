
#' apply tmle_single_p for num_cf times
#'
#' @param data similar as aipw_single_p() function
#' @param exposure similar as aipw_single_p() function
#' @param outcome similar as aipw_single_p() function
#' @param covarsT similar as aipw_single_p() function
#' @param covarsO similar as aipw_single_p() function
#' @param learners similar as aipw_single_p() function
#' @param control similar as aipw_single_p() function
#' @param num_cf number of partitions
#' @param n_split similar as aipw_single_p() function
#' @param rand_split logical value; if be TRUE, discordant splits for exposure and outcome model are chosen at random ; otherwise chosen systematically.
#' @param seed numeric value to reproduce the splits
#' @return a tibble of the estimates
#' @export
#'
#' @examples
#' sum(1:5)
#'
#'
tmle_multiple_p <- function(data, exposure, outcome, covarsT, covarsO, learners, control, num_cf, n_split, rand_split = TRUE, seed = 145){

  #Initialize results
  runs <- tibble(r1=double(), r0=double(), rd=double(), v1=double(), v0=double(), vd=double())

  #Run on num_cf splits
  for(cf in 1:num_cf){
    runs <- bind_rows(runs, DC_TMLE_Single_p(data, exposure, outcome, covarsT, covarsO, learners, control, n_split))

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

