\name{shypFun.02110}
\alias{shypFun.02110}
\title{Kosugi-Mualem Model (2 Parameter Model)
}
\description{Calculates the soil hydraulic property function values based on given pressure heads
}
\usage{ shypFun.02110(p, h)
}
\arguments{
\item{p}{vector of the 6 Kosugi-Mualem model parameters, order is sensitve and has to be given as:
\tabular{lll}{
\code{thr}\tab{residual water water content [cm cm-3]}\cr
\code{ths}\tab{saturated water water content [cm cm-3]}\cr
\code{hm}\tab{air entry pressure head [cm]}\cr
\code{sigma}\tab{width of pore size distribution [ - ]}\cr
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
\bold{ Kosugi, K.}: Lognormal distribution model for unsaturated hydraulic properties, 
Water Reourc. Res., 32(9), 2697-2703, <doi:10.1029/96WR01776>, 1996.
}
\author{
Tobias KD Weber
}

\examples{
p <- c("thr" = 0.1, "ths" = 0.4, "hm" = 100, "sigma" = 2, "Ks" = 100, "tau" = .5)
h <- 10^seq(-2, 6.8, length = 197)
shyp.L <- shypFun.02110(p, h)
}
\keyword{ ~Kosugi}
