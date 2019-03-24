\name{shypFun.01110}
\alias{shypFun.01110}
\title{ van Genuchten-Mualem Function
}
\description{Calculates the soil hydraulic property function values based on given pressure heads
}
\usage{shypFun.01110(p, h)
}
\arguments{
\item{p}{vector of the 6 van Genuchten-Mualem model parameters, order is sensitve and has to be given as:
\tabular{lll}{
\code{thr}\tab{residual water water content [cm cm-3]}\cr
\code{ths}\tab{saturated water water content [cm cm-3]}\cr
\code{alf1}\tab{van Genuchten alpha [cm-3]}\cr
\code{n1}\tab{van Genuchten n [-]}\cr
\code{Ks}\tab{saturated conductivity [cm d-1]}\cr
\code{tau}\tab{exponent of \code{Se} in the capillary conductivity model, sometimes denoted in the literature as \code{l} [-]}

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
\bold{van Genuchten, M.T.}: Closed-form equation for predicting the hydraulic conductivity of 
 unsaturated soils, Soil Sci Soc Am J, 44(5), 892-898, <doi:10.2136/sssaj1980.03615995004400050002x>, 1980.
}
\author{Tobias KD Weber
}

\examples{
p <- c(0.1, 0.4, 0.01, 2, 100, .5)
h <- 10^seq(-2, 6.8, length = 197)
shyp.L <- shypFun.01110(p, h)
}

\keyword{ ~VGM}
