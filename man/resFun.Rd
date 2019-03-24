\name{resFun}
\alias{resFun}
\title{Calculation of the Objective Function Value
}
\description{Contains the objective functions to calculate (weighted) sum of squared residuals or likelihoods.
The assumption made is that the residuals are identically, independantly and normally distributed.
The normal distribution of the model residuals is standardly given with mean = 0, and variance = 1.if weighting is not considered.
There are three methods to consider weights: a) fixed skalar values for each data type, b) a vector of weights for each 
data type. The vector has to have the same length as the vector of the data type.
c) It is possible to estimate the 

}
\usage{
resFun(p, shpmodel = "01110", retdata = NULL, condata = NULL, pretrans = NULL,
      weight = NULL, method = "rss", 
      trim.query = FALSE, ivap.query = NULL, hclip.query = FALSE, parL = NULL)

}
\arguments{
\item{p}{Vector of model parameters handed used to calculate the soil hydraulic property model values in \link[=shypFun]{shypFun}.
Depends on \code{shpmodel} and the pressure head values specified in \code{retdata} and \code{condata} 
}
\item{shpmodel}{Character identifying the soil hydraulic property model. See \link[=shypFun]{shypFun}.
}
\item{retdata}{a dataframe or matrix with 2 columns. The first with pressure head values in [cm] and the second with volumetric water contents in [cm cm-3].
}
\item{condata}{a dataframe or matrix with 2 columns. The first with pressure head values in [cm] and the second with hydraulic conductivity values in log10[cm d-1].
}
\item{pretrans}{a vector to back transform the parameters before the soil hydraulic property function values calculated.
}
\item{weight}{specification of weight method. See \link[=weightFun]{shypFun}

}
\item{method}{
 \tabular{lll}{
"rss" default for the optimisation algortihm DEoptim. resFun returns skalar sum of squared (weighted) residuals. \cr
"res" resFun returns a list with vectors of weighted residuals. Required for post hoc analyses. \cr
}
}
\item{trim.query}{ default FALSE. If a trimodal soil hydraulic property function is used, this has to be specified by setting the argument to (\emph{TRUE}) which ensures the sum of modal weights == 1.
}
\item{ivap.query}{Default is FALSE, and no\emph{ivap} method is specified. See \link[=KvapFun]{KvapFun}.
}
\item{hclip.query}{Implemented purely for future compatability. Currently no use.
}
\item{parL}{defaults to \code{NULL}, only inserted for compatbility with \code{modMCMC} used in \link[=shypEstFun]{shypEstFun}. \code{modMCMC}, only
handled parameters which are estimated, other model parameters need to be passed through parL. See Details of \link[=shypEstFun]{shypEstFun}.
}
}
\details{ Model errors may be specified or estimated as nuisance parameters weighting the data classes. In case the model error !=1, the output statistics are weighted accordingly. 

 \tabular{lll}{
\code{user}\tab{user defined weights}\cr
\code{none}\tab{no weights are considered, i.e. no measurement error assumed}\cr
\code{range}\tab{rescaling (normalization of observations to the intervall [0,1]}\cr
\code{fix1}\tab{fixed scalar weight for THETA is 1/0.05^2 and weight for log10K is 1 }\cr
\code{est1}\tab{Two scalar model weights 1/sigma_i^2 are treated as free parameters to be estimated by inversion, one for \code{THETA} and one for \code{log10K}. Only simultaneously estimateable}\cr
} 
 }
\value{Returns scalar of sum of squared (weighted) residuals or vector of weighted residuals, as specified by 

 \tabular{lll}{
\code{rss}\tab{scalar sum of squared (weighted) residuals}\cr
\code{res}\tab{vector of weighted residuals}\cr
\code{loglik}\tab{log-likelihood value}\cr
} 

}
\references{Example and details on weighted objective functions using multimodal soil hydraulic property functions
{{\bold{Weber, T.K.D., Iden, S.C., and W. Durner} Unsaturated hydraulic properties of Sphagnum moss and peat reveal trimodal pore-size distributions, Water Resour Res, <doi:10.1002/2016WR019707>, 2017.}}
}
\author{
Tobias KD Weber
}
\note{ Calculates the objective function value as the sum of squared weighted model errors.

}



\examples{
# load data
data("shpdata1")

# observations
retdata <- shpdata1$LFH1$wrc[!is.na(shpdata1$LFH1$wrc[,1]),]
condata <- shpdata1$LFH1$hcc

# 7 - resFun ------------------------------------------------------------
# soil hydraulic property model parameters, van Genuchten-Mualem
p <- c("thr" = 0.16, "ths" = 0.46, "alf1" = 0.03, "n1" = 1.42, "Ks" = 26, "tau" = .5)

# calculate weighted residuals
wres <- resFun(p, retdata = retdata, condata = condata, pretrans = NULL,
weight = list("wth" = 0.0025, "wKh" = 1), method = "res", trim = FALSE)

## residuals of the soil water retention curve [-]
theta.wres <- wres[1:dim(retdata)[1]]

## residuals of the log10 hydraulic conductivity curve [cm/d]
log10K.wres <- wres[(dim(retdata)[1]+1) : length(wres)]

}
\keyword{ ~resFun}
