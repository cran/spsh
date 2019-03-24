\name{sncFun}
\alias{sncFun}
\title{Non-capillary Saturation Function
}
\description{The general purpose method to calculate the effective non-capillary saturation is directly obtained from any arbritrary expression for the 
rescaled capillary saturation function as described by Weber et al. (2019). Examples of capillary saturation functions are the well known van Genuchten (1980), Fredlund and Xing (1993), and Kosugi (1996) functions.
}
\usage{sncFun(h, scap)
}
\arguments{
  \item{h}{A vector of \code{n} pressure head values for which \code{scap} was calculated}
  \item{scap}{vector of \code{n} monotonically decreasing capillary saturation function values calculated by \code{shypFun}, rescaled between 0 and 1.
   }
  
}
\details{
More details in Weber et al. (2019)
}
\value{ 
    \item{snc}{A vector of \code{n} element with calculated volumetric moisture content}
}
\references{
Fredlund D.G., and Xing, A.: Equations for the soil-water characteristic curve. Can. Geotech. J. 31:521-532, 
<doi:10.1139/t94-061>, 1994.

Kosugi, K. Lognormal distribution model for unsaturated hydraulic properties, 
Water Reourc. Res., 32(9), 2697-2703, <doi:10.1029/96WR01776>, 1996.

van Genuchten, M.T.: Closed-form equation for predicting the hydraulic conductivity of 
unsaturated soils, Soil Sci Soc Am J, 44(5), 892-898, <doi:10.2136/sssaj1980.03615995004400050002x>, 1980.
            
Weber, T.K.D., Durner, W., Streck, T., and Diamantopoulos, E.: A modeluar framework for modelling unsaturated soil hydraulic properties over the full moisture range, in revision, 2019.
}
\author{
Tobias KD Weber
}
\note{The function requires a numerical solution to the integral of Eq. 6 in Weber et al. (2018) and; therefore, it is advisable to use a sufficient 
amount of pressure head data values to minimise the numerical error.
}

\examples{
# set variables
p <- c(0.1, 0.4, 0.01, 2, 100, .5)
h <- 10^seq(-2, 6.8, length = 197)

# Calculate the capillary and non-capillary saturation function.
Se <- shypFun(p, h, shpmodel = "01110")$Se     
Snc <- sncFun(Se)
  }
\keyword{ ~VGM}
