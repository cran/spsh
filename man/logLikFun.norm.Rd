\name{logLikFun.norm}
\alias{logLikFun.norm}
\title{Calculation of the Log-likelihood assuming Identially, independenzly and Normally Distributed errors
}
\description{Calculates the i-th log-likelihood of each \emph{y-yhat} pair as described in Seber and Wild (2003, ISBN:9780471617600) Eq. 2.46.
}
\usage{logLikFun.norm(y, yhat, sigma)
}
\arguments{
\item{y}{A vector of \code{n} observed properties/variables of interest.
}
\item{yhat}{A vector of \code{n} model simulated properties/variables of interest.
}
\item{sigma}{A vector of length 1 considering homoscedastic residuals.
}

}
\details{
The underlying assumption is, that the model residuals (errors) are independently, and identically distributed (i.i.d.) following a normal distribution.
}

\value{ \emph{log-likelihood} value of an normal distribution with N~(0, \emph{sigma^2})

}

\references{
\bold{Seber, G.A.F. and Wild, C.J.}: Nonlinear Regression, Wiley-Interscience, New Jersey, 2003.
}
\author{
Tobias K.D. Weber <tobias.weber@uni-hohenheim.de>,
}
\note{
The assumption of i.i.d. and normal distribution is best investigated \emph{a posteriori}.
}
\examples{

# homoscedastic residuals
sig.s <- .01
y.scat <- rnorm(100, 0, sig.s)

yhat <- (1:100)^1.2
y <- yhat+y.scat
sum(logLikFun.norm(y, yhat, sig.s))
plot(yhat-y)


}
\keyword{ ~PTF}
 