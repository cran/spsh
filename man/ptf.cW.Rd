\name{ptf.cW}
\alias{ptf.cW}
\title{Corrected Weynants et al. (2009) Pedotransfer Function
}
\description{This function predicts van Genuchten-Mualem model parameters, setting the residual water
content to zero.
}
\usage{ptf.cW(CLAY, SAND, BD, OC)
}
\arguments{
\item{CLAY}{A vector of \code{n} elements with soil clay content (particle diameters <= 2 x 10e-6 m), in percent [0, 100]).
}
\item{SAND}{A vector of \code{n} elements with soil sand content (particle diameters < 2 mm and > 50 x 10e-6 m), in percent [0, 100]).
}
\item{BD}{A vector of \code{n} elements with soil bulk density (g/cm3).
}
\item{OC}{A vector of \code{n} elements with soil organic carbon content, in percent [0, 100].
}
}
\details{Pedotransfer function returns the van Genuchten - Mualem model parameters given \code{CLAY}, \code{SAND}, \code{BD}, and \code{OC}.
The correction of the original paper presented by Weynants et al (2009), were made by Weihermueller et al., (2017), which is implemented.
}

\value{ 
List with the van Genuchten-Mualem parameters, each as a vector of \code{n} elements
\tabular{lll}{
\code{thr}\tab{Residual water content (-), alway equal to zero}\cr
\code{ths}\tab{Saturated water content (-)}\cr
\code{alf1}\tab{Shape parameter (cm^-1)}\cr
\code{n1}\tab{Shape parameter (-)}\cr
\code{K0}\tab{Hydraulic conductivity at 0 potential (cm/day)}\cr
\code{tau}\tab{Shape parameter (-)}\cr
}
}

\references{
\bold{Weynants, M., Vereecken, H., and Javaux, M. M.}: Revisiting Vereecken pedotransfer functions: Introducing a closed-form hydraulic model.
Vadose Zone J, 8(1), 86-95, <doi:10.2136/vzj2008.0062>, 2009.

\bold{Weihermueller, L., Herbst, M., Javaux, M., and Weynants, M.}: Erratum to "Revisiting Vereecken Pedotransfer Functions: Introducing a Closed-Form Hydraulic Model", Vadose Zone J, 16(1), <doi:10.2136/vzj2008.0062er>, 2017. 
}
\author{
Melanie Weynants <mweynants@gmail.com>
Tobias K.D. Weber <tobias.weber@uni-hohenheim.de>
}
\note{
The PTF is not suitable for predicting thehydraulic conductivity curve at pressured heads > -6 cm. (Weynant et al., 2009)
}
\examples{
ptf.cW(CLAY = .4, SAND = .4, BD = 1.6, OC = .5)
}
\keyword{ ~PTF}
 