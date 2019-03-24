\name{shypFun.01210}
\alias{shypFun.01210}
\title{Durner Model (bimodal van Genuchten Mualem)
}
\description{Calculates the soil hydraulic property function values based on given pressure heads
}
\usage{shypFun.01210(p, h)
}
\arguments{
\item{p}{vector of the 9 bimodal van Genuchten-Mualem model parameters, order is sensitve and has to be given as:
\tabular{lll}{
\code{thr}\tab{residual water water content [cm cm-3]}\cr
\code{ths}\tab{saturated water water content [cm cm-3]}\cr
\code{alf1}\tab{van Genuchten alpha [cm-3]}\cr
\code{n1}\tab{van Genuchten n [-]}\cr
\code{w1}\tab{fraction of the first modality [-], \code{w2} is internally computed as \code{w2 = 1-w1}}\cr
\code{alf2}\tab{van Genuchten alpha of the second modality[cm-3]}\cr
\code{n2}\tab{van Genuchten n of the second modality [-]}\cr
\code{Ks}\tab{saturated conductivity [cm d-1]}\cr
\code{tau}\tab{exponent of \code{Se} in the capillary conductivity model, sometimes denoted in the literature as \code{l} [-]}\cr

}}
\item{h}{pressure heads [cm] for which the corresponding retention and conductivity values are calculated
}
}
\details{The function solves analytically the spec. water capacity function and integral to the capillary bundle model 
}
\value{ returns a \code{list} with calculations at specified \code{h}:
\item{theta}{calculated volumetric moisture content}
\item{Se}{calculated saturation}
\item{cap}{specific water capacity function}
\item{psd}{pore size distribution}
\item{Kh}{Hydraulic conductivity values}
}
\references{
 {\bold{Durner, W.}: Hydraulic conductivity estimation for soils with heterogeneous pore structurre, WRR, 30(2), 211-223, <doi:10.1029/93WR02676>, 1994.}
 
 {\bold{Priesack, E. and Durner, W.}: Closed-Form Expression for the Multi-Modal Unsaturated Conductivity Function. Vadose Zone Journal 5:121-12, <doi:10.2136/vzj2005.0066>, 2006.}\cr

}
\author{
Tobias KD Weber
}

\examples{
p <- c("thr" = 0.1, "ths" = 0.4, alf1 = 0.5, "n1" = 3,
"w1" = .6, "alf2" = 0.01, "n2" = 1.6, 
"Ks" = 100, "tau" = .5)
h <- 10^seq(-2, 6.8, length = 197)
shyp.L <- shypFun.01210(p, h)
}

\keyword{ ~VGM}
