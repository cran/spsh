\name{shypFun.03110}
\alias{shypFun.03110}
\title{Unimodal Fredlund-Xing Model
}
\description{Calculates the soil hydraulic property function values based on given pressure heads. The function calculates the base function of Fredlund and Xing.

}
\usage{shypFun.03110(p, h)
}
\arguments{
\item{p}{vector of the 6 Fredlund-Xing model parameters, order is sensitve and has to be given as:
\tabular{lll}{
\code{thr}\tab{residual water water content [cm cm-3]}\cr
\code{ths}\tab{saturated water water content [cm cm-3]}\cr
\code{alf1}\tab{inverse of the air entry pressure head [cm]}\cr
\code{n1}\tab{width of pore size distribution [ - ]}\cr
\code{Ks}\tab{saturated conductivity [cm d-1]}\cr
\code{tau}\tab{exponent of \code{Se} in the capillary conductivity model, sometimes denoted in the literature as \code{l} [ - ]}\cr

}}
\item{h}{pressure heads [cm] for which the corresponding retention and conductivity values are calculated
}
}

\details{The function numerically solves the spec. water capacity function and itnegral to Mualem's conductivity model.


}
\value{ returns a \code{list} with calculations at specified \code{h}:
\item{theta}{calculated volumetric moisture content}
\item{Se}{calculated saturation}
\item{cap}{specific water capacity function}
\item{psd}{pore size distribution}
}
\references{
\bold{Fredlund D.G., and A. Xing.}: Equations for the soil-water characteristic curve, Can. Geotech. J., 31:521-532, <doi: 10.1139/t94-061>,1994.
}
\author{
Tobias KD Weber
}

\examples{
p <- c(0.1, 0.4, 0.01, 2, 100, .5)
h <- 10^seq(-2, 6.8, length = 197)
shyp.L <- shypFun.03110(p, h)
}

\keyword{ ~Fredlund}
