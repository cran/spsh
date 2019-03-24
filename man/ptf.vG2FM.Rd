\name{ptf.vG2FM}
\alias{ptf.vG2FM}
\title{Parameter Tranfer Function for Weber et al. (2019) model
}
\description{Predicts Weber et al. (2019) model parameters in the van Genuchten-Mualem variant \emph{01110FM}, from previously obtained van Genuchten-Mualem parameters for the constrained van Genuchten-Mualem model.
}
\usage{ptf.vG2FM(x)
}
\arguments{
\item{x}{A vector of 6 van Genuchten-Model parameters.
}
}
\details{Pedotransfer function returns the van Genuchten - Mualem model \emph{01110} parameters in the Brunswick-Model variant \emph{01110FM}, based on previously determined van Genuchten-Mualem parameters. The transfer function is based on an ordinary linear regression between the i-th \emph{01110} and \emph{01110FM}. The paraemeters
were based on model fits to a dataset of ~200 samples with retention and conductivity measurements.}

\value{ 
vector of van Genuchten-Mualem model parameters, the order of which is sensitive.
\tabular{lll}{
\code{thr}\tab{Residual water content (-), alway equal to zero}\cr
\code{ths}\tab{Saturated water content (-)}\cr
\code{alf1}\tab{Shape parameter (cm^-1)}\cr
\code{n1}\tab{Shape parameter (-)}\cr
\code{Ks}\tab{Hydraulic conductivity at 0 potential (cm/day)}\cr
\code{tau}\tab{Shape parameter (-)}\cr
}
}

\references{
\bold{Weber, T.K.D, and Diamantopoulos, E.}: A simple model to derive the framework soil hydraulic property model parameters from known van Genuchten-Mualem model parameters, unpublished, 2018.}
\author{
Efstathios Diamantopoulos <ed@plen.ku.dk> 
Tobias K.D. Weber <tobias.weber@uni-hohenheim.de>
}
\note{ The parameter transfer function was derived by ordinary linear regression, by regressing the estimated the Weber et al. (2018) framework model parameters Genuchten Mualem variant frommeasured soil water retention and hydraulic conductivity data against the constrained van Genuchten model parameters. 
The regression of \code{alf1} and \code{Ks}, and \code{(n-1)} was done in the \code{log[10])}-transformed space,
and \code{Kncs} is predicted as the fraction w of Ks
}
\examples{
p = c(0.08, 0.42, 0.01, 1.5, 100, 0.5)
ptf.vG2FM(p)
}
\keyword{ ~PTF}
 