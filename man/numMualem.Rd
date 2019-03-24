\name{numMualem}
\alias{numMualem}
\title{Function to Numerically Compute the Mualem Integral
}
\description{This function will calculate Mualems integral and returns hydraulic conductivity values.

}
\usage{
numMualem(h, pcon, scap)
}

\arguments{
\item{h}{vector of length \code{l} with pressure head values.
}
\item{pcon}{vector of soil hydraulic conductivitymodel parameters, the first argument has to be the tortuosity parameter 'tau'

}
\item{scap}{vector of length \code{l} of calculated (capillary) saturation values 
}
}

\details{
The numerical solution of Mualems integral relies on the trapezoidal rul of integration.
 }

\value{
\item{numMualem}{returns a vector of length \code{l} with calulcated conductivity values.
}
}
\references{
\bold{Mualem, Y.}:A New Model for Predicting the Hydraulic Conductivity of Unsaturated Porous Media,Water Resources Research, 12(3), <doi:10.1029/WR012i003p00513>, 1976.
}

\author{ Tobias KD Weber

}

\examples{

h <- 10^seq(-3, 6.8, length = 501)
p = c(.05, .5, .01, 1.8, 100, .5)
shyp.L <- shypFun.01110(p, h)

Ks <- p[5]
tau <- p[6]
Se <- shyp.L[['Se']]
Khrnum <- numMualem(h, pcon = tau, scap = Se) 

Khnum <- Ks * Se^tau * Khrnum

plot(log10(h), log10(shyp.L[['Kh']]), ylim = c(-10, 2.3), 
xlim = c(-1,6), ylab = "log10 Kunsat [ cm/d ]", xlab = "pF [ - ]", type = "l", lwd = 8)
lines(log10(h), log10(Khnum), col = "red", lwd = 2)

 
}

\keyword{ ~Kvap }
\keyword{ ~kvap }
