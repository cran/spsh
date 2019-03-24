\name{weightFun}
\alias{weightFun}
\title{Specification of Weights for the Data Groups Retention Data and Conductivity Data.
}
\description{Weights can be fixed to suggested standards, fixed by the user, or estimated as 
additional nuisance parameters.

}
\usage{
weightFun(weightmethod = "fix1", retdata, condata, parL = NA)
}

\arguments{
 \item{weightmethod}{character specifying the method of selecting model weights
 }
 \item{retdata}{a dataframe or matrix with 2 columns. The first with pressure head values in [cm] and the second with volumetric water contents in [cm cm-3].
 }
 \item{condata}{a dataframe or matrix with 2 columns. The first with pressure head values in [cm] and the second with hydraulic conductivity values log10[cm d-1].
 }
 \item{parL}{ Defaults to \code{NA} has to be provided if \code{weightmethod == "est1"} . See Details of (\link[=shypEstFun]{shypEstFun} for explanation of \code{parL})} 
}

\details{Character specifying \code{weightmethod}

\tabular{lll}{
\code{user}\tab{user defined weights}\cr
\code{none}\tab{no weights are considered, i.e. no measurement error assumed}\cr
\code{range}\tab{rescaling (normalization of observations to the intervall [0,1]}\cr
\code{fix1}\tab{fixed scalar weight for THETA is 0.05^2 and weight for log10K is 1 }\cr
\code{fix2}\tab{vector with the length of number of observations as given in \code{retdata} and \code{condata} are given, fixed to weight for THETA is 0.05^2 and weight for log10K is 1 }\cr
\code{est1}\tab{Two scalar model weights (sigma^-2) are treated as free parameters to be estimated by inversion, one for \code{THETA} and one for \code{log10K}}\cr
}
}

\value{The function returns a \code{list} of \code{weights} as specified through \code{weightmethod}

If \code{weightmethod} is set to \code{est1} and \code{parL} is given as an extra argument, the function returns a list wich is concatenated to the \code{parL} used in \code{shypEstFun}
providing extra information on the nuisance parameters. Alternatively, \code{parL} can be passed as an argument to \code{shypEstFun} directly, accounting for the two additional nuisance parameters at the end of the respective vectors.
}
\author{
Tobias KD Weber
}

\examples{
# Example 1 | fixed weights
weight.fix.L <- weightFun("fix1") 

## Example 2 | range of measure data
data(shpdata1)

wrc <- shpdata1$TS1$wrc
hcc <- shpdata1$TS1$hcc
# Remove NAs
hcc <- shpdata1$TS1$hcc[!is.na(shpdata1$TS1$hcc[,1] ),]
weight.fix.L <- weightFun("range", wrc, hcc) 

}

\keyword{ ~weights}
