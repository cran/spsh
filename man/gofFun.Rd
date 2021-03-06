% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gofFun.r
\name{gofFun}
\alias{gofFun}
\title{Goodness-of-fit and Information Criteria}
\usage{
gofFun(
  phat,
  shpmodel = "01110",
  retdata = NULL,
  condata = NULL,
  weight,
  psel,
  ivap.query = NULL,
  hclip.query = FALSE
)
}
\arguments{
\item{phat}{Vector of non-transformed (back-transformed) model parameters after estimation, i.e. the best fit or maximum likelihood estimate.}

\item{shpmodel}{Character specifying the soil hydraulic property model.}

\item{retdata}{Dataframe or matrix with 2 columns. The first with pressure head values in \ifelse{html}{\out{log<sub>10</sub>}}{\eqn{log_10}} [cm], i.e. pF values, and the second with volumetric water contents in [cm cm-3].}

\item{condata}{Dataframe or matrix with 2 columns. The first with pressure head values in \ifelse{html}{\out{log<sub>10</sub>}}{\eqn{log_10}} [cm], i.e. pF values, and the second with hydraulic conductivity values \ifelse{html}{\out{log<sub>10</sub>}}{\eqn{log_10}} [cm d-1].}

\item{weight}{List of the model residual weights used or estimated by the parameter estimation scheme, to calculate the weighted statistical analyses.}

\item{psel}{Vector specifying the selected parameters for the parameter estimation from \code{parL}.}

\item{ivap.query}{Specification of \emph{ivap} method, if FALSE selected, no isothermal vapour conductivity is consideredSee \code{Details}}

\item{hclip.query}{Implemented purely for future compatability. Currently no use. See \code{Details}}
}
\description{
Calculates goodness-of-fit criteria and the likelihood-based Akaike 
and Bayesian Information Criterion based on a given parameter set,
typically from the estimation procedure.
}
\details{
Output for data groups.
 \tabular{llll}{
    \code{th}\tab      {\code{list} with goodness of fit statistics for the retention curve \emph{see below}}\cr
    \code{logKh}\tab   {\code{list}with output same as \code{th} but for the \ifelse{html}{\out{log<sub>10</sub>}}{\eqn{log_10}} fitted conductivity curve}\cr
    \code{combined}\tab{\code{list} with \code{AIC}, \code{AICc}, and \code{BIC} calculated for the multi-objective function if arguments \code{retdata} and \code{condata} are both \code{!NULL}}\cr
   }
Statistical analyses of the inverse modelling results.
  \tabular{lll}{
    \code{me}\tab{mean (weighted) error}\cr
    \code{mae}\tab{mean absolute (weighted) error}\cr
    \code{mse}\tab{mean squared (weighted) error}\cr
    \code{rss}\tab{sum of squared (weighted) errors}\cr
    \code{rmse}\tab{root mean squared (weighted) error}\cr
    \code{AIC}\tab{Akaike Information Criteria}\cr
    \code{AICc}\tab{corrected Akaike Information Criteria}\cr
    \code{BIC}\tab{Bayesian Information Criteria}\cr
    \code{m}\tab{number of observations}\cr
            }
}
\examples{
data("shpdata1")
retdata <- shpdata1$TS1$wrc
condata <- shpdata1$TS1$hcc
condata <- condata[!is.na(condata[,1]),]
# Parameter list
parL <- list("p" = c("thr"= 0.05, "ths" = 0.45, "alf1" = 0.01, "n" = 2, "Ks" = 100, "tau" = .5),
             "psel" = c(1, 1, 0, 1, 1, 1),
             "plo" = c(0.001 , 0.2, 0.001, 1.1, 1, -2),
             "pup" = c(0.3, 0.95, 1, 10, 1e4, 10)
             )
# Calulation of the goodness of fit.
gofL <-gofFun(parL$p, shpmodel = "01110", retdata = retdata, condata = condata, 
              weight = weightFun(weightmethod = "fix1"), parL$psel, 
              ivap.query = NULL, hclip.query = FALSE)
}
\references{
\insertRef{Hoege.2018}{spsh}
}
\author{
Tobias KD Weber , \email{tobias.weber@uni-hohenheim.de}
}
