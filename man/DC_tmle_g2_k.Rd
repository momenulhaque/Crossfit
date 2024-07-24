% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DC_tmle_g2_k.R
\name{DC_tmle_g2_k}
\alias{DC_tmle_g2_k}
\title{Estimate Average Treatment Effect (ATE) using from TMLE estimator using cross-fit algorithm (generalization 2)}
\usage{
DC_tmle_g2_k(
  data,
  exposure,
  outcome,
  covarsT,
  covarsO,
  family.y = "binomial",
  learners,
  control,
  n_split = 3,
  num_cf = 10,
  rand_split = FALSE,
  gbound = 0.025,
  alpha = 1e-17,
  seed = 146,
  conf.level = 0.95
)
}
\arguments{
\item{data}{a data frame of tibble}

\item{exposure}{name of exposure variable}

\item{outcome}{name of outcome variable}

\item{covarsT}{a vector of names of covaraites for treatment model}

\item{covarsO}{a vector of names of covaraites for outcome model}

\item{family.y}{it is the family for outcome model. It can \code{binomial() (default)} or \code{"gaussian"}}

\item{learners}{similar as\code{SL.library()} in \code{SuperLearner} package.}

\item{control}{similar as \code{cvControl()} in \code{SuperLearner} package.}

\item{n_split}{number of splits used, default \code{n_split = 3}}

\item{num_cf}{number of repetition done. The default is 5.}

\item{rand_split}{logical value; if be FALSE \code{(default)}, discordant splits for exposure and outcome model are chosen systematically; otherwise chosen randomly.}

\item{gbound}{value between (0,1) for truncation of predicted probabilities. The defaults are 0.025 and 0.975. See \code{tmle::tmle()} for more information.}

\item{alpha}{used to keep predicted initial values bounded away from (0,1) for logistic fluctuation. The defaults are 1e-17 and 1-1e-17.}

\item{seed}{numeric value to reproduce the splits distribution}

\item{conf.level}{confidence limit for confidence interval, \code{default = 0.95}.}
}
\value{
It return a list of two elements. The first element \code{ATE} is a tibble of the estimates. The \code{weight} is a tibble of weights of learners for two different models.
}
\description{
Estimate Average Treatment Effect (ATE) using from TMLE estimator using cross-fit algorithm (generalization 2)
}
\examples{

# See the README file for details

sum(1:5)


}