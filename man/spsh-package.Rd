\name{spsh-package}
\alias{spsh-package}
\alias{spsh}
\docType{package}
\title{Soil Physics and Soil Hydrology (spsh) Package with Options for Parameter Estimation 
}
\description{
This package provides a multitude of spsh models, in particular the soil hydraulic property (SHP) models; i.e. the soil water retention curve (WRC) and unsaturated hydraulic conductivity curve (HCC). For the WRC and HCC, (weighted) parameter estimation options exist. Global optimisation is done with the Differential Evolution algorithm \emph{DEoptim}. A number of goodness-of-fit metrics are implemented, some of which routed in Bayesian theory.
Also the \emph{DRAM} Markov-chain Monte Carlo algorithm as implemented in \code{modMCMC}, can be used to to sample the posterior. 
A pedotransfer function to predict soil hydraulic property model parameters based on soil textural information is included.

The analyses of the goodness-of-fit is based on the best parameter estimate. Note, that non-uniform model priors can currently not be specified.

}
\details{
\tabular{ll}{
Package: \tab spsh\cr
Type: \tab Package\cr
Version: \tab 1.0.4\cr
Date: \tab 2019-03-19\cr
License: \tab GPL-2\cr
}
\bold{Supported Soil Hydraulic Property Models} 
 \tabular{lll}{
 \code{01110}\tab\tab{unimodel van Genuchten-Mualem model, with the contraint of m = 1-1/n (van Genuchten, 1980)}\cr
 \code{01210}\tab\tab{bimodel van Genuchten-Mualem model, with the contraint of \eqn{m_i = 1-1/n_i}(Durner, 1994)}\cr
 \code{01310}\tab\tab{trimodal van Genuchten-Mualem model, with the contraint of \eqn{m_i = 1-1/n_i} (Durner, 1994)}\cr
 \code{02110}\tab\tab{unimodal 2 parametric Kosugi-Mualem model (Kosugi, 1996)}\cr
 \code{03110}\tab\tab{unimodal base Fredlund-Xing-Mualem model, with the contraint of m = 1-1/n (Fredlund and Xing, 1994)}\cr
}
\bold{Framework Model (Brunswick Model)} 

The Framework Model, as proposed by Weber et al. (2019), can be used with the implemented soil hydraulic property functions.

This option can be accessed by adding "FM" to the code specified in \emph{shpmodel}, e.g. in \link[=shypFun]{shypFun}


\bold{Supported Pedo-transfer function } 
 \tabular{lll}{
\code{cW}\tab\tab{ corrected Weynant model of Weynants et al. (2009) and Weihermueller et al. (2017)}
}
 
(For a mathematical description refer to original literature in the references or the provided R code)
}
\author{
Tobias K.D. Weber\email{tobias.weber@uni-hohenheim.de} \cr
 \cr
Biogeophysics Department \cr
Institute of Soil Science and Land Evaluation\cr
University of Hohenheim \cr
Germany

}
\references{
\bold{Durner, W.}: Hydraulic conductivity estimation for soils with heterogeneous pore structurre, WRR, 30(2), 211-223, <doi:10.1029/93WR02676>, 1994.\cr

\bold{Fredlund D.G., and Xing, A.}: Equations for the soil-water characteristic curve, Can. Geotech. J., 31:521-532, <doi:10.1139/t94-061>,1994.\cr

\bold{Kosugi, K.}: Lognormal distribution model for unsaturated hydraulic properties, Water Reourc. Res., 32(9), 2697-2703, <doi:10.1029/96WR01776>, 1996.\cr

\bold{Mualem, Y.}: New model for predicting hydraulic conductivity of unsaturated porous media. Water Resour. Res. 12(3): pp. 513-522. <doi:10.1029/WR012i003p00513>, 1976\cr

\bold{Priesack, E. and Durner, W.}: Closed-Form Expression for the Multi-Modal Unsaturated Conductivity Function, Vadose Zone Journal 5:121-12, <doi:10.2136/vzj2005.0066>, 2006.\cr

\bold{van Genuchten, M.T.}: Closed-form equation for predicting the hydraulic conductivity of unsaturated soils. Soil Sci Soc Am J 44(5): pp. 892-898. <doi:10.2136/sssaj1980.03615995004400050002x>, 1980\cr

\bold{Weber, T.K.D., Durner, W., Streck, T., and Diamantopoulos, E.}: A modular framework for modelling unsaturated soil hydraulic properties over the full moisture range, WRR , in rev, 2019.\cr

\bold{Weihermueller, L., Herbst, M., Javaux, M., and Weynants, M.}: Erratum to "Revisiting Vereecken Pedotransfer Functions: Introducing a Closed-Form Hydraulic Model", Vadose Zone J, 16(1), <doi:10.2136/vzj2008.0062er>, 2017.\cr

\bold{Weynants, M., Vereecken, H., and Javaux, M. M.}: Revisiting Vereecken pedotransfer functions: Introducing a closed-form hydraulic model. Vadose Zone J, 8(1), 86-95, <doi:10.2136/vzj2008.0062>, 2009.\cr
}
\seealso{
\link[=shypFun]{shypFun}, \link[=shypEstFun]{shypEstFun}
}

\examples{
p <- c("thr" = 0.1, "ths" = 0.4, "alf1" = 0.01, "n1" = 2, "Ks" = 100, "tau" = .5)
h <- 10^seq(-2, 6.8, length = 197)
shyp.L <- shypFun.01110(p, h)
}

\keyword{ ~soil hydrology }
\keyword{ ~soil physics }
\keyword{ ~parameter estimation}
