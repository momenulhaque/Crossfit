#' A simulated data set having 3000 rows and 7 columns.
#'
#' The `statin_sim_data` is generated from the guidelines described by `Zivich and Breskin (2021)`
#'
#'
#' @format ## `statin_sim_data`
#' A data frame with 3000 rows and 7 following columns:
#' \describe{
#'   \item{Y}{Atherosclerotic cardiovascular disease (ASCVD) status}
#'   \item{statin}{exposure of statin usage}
#'   \item{age}{patient's age}
#'   \item{ldl_log}{log of low-density lipoprotein measurement}
#'   \item{diabetes}{Diabetes status}
#'   \item{risk_score}{Risk score of ASCVD for each patient}
#'   \item{risk_score_cat}{Risk category of ASCVD for each patient. The patient belongs to 3 category has the greatest risk}
#'
#' }
#'
#'
#' @references 1. Paul N Zivich and Alexander Breskin (2021),
#' Machine learning for causal inference: on the use of cross-fit estimators.
#' \emph{Epidemiology}, 32(3):393–401.
#' \doi{10.1097/EDE.0000000000001332}
#'
#'
"statin_sim_data"

