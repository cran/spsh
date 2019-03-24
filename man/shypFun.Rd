\name{shypFun}
\alias{shypFun}
\title{Wrapper Function for all Supported Soil Hydraulic Property Models.
}
\description{ This function allows to select a soil hydraulic proerty model.
}
\usage{shypFun(p ,h , shpmodel = "01110", ivap.query = NULL)
}

\arguments{
\item{p}{Vector of the 6 van Genuchten-Mualem model parameters, order is sensitve. cf respective model documentation}
\item{h}{Pressure heads [cm] for which the corresponding retention and conductivity values are calculated
}
\item{shpmodel}{ character
 \tabular{lll}{
\code{01110}\tab{unimodel van Genuchten-Mualem model, with the contraint of m = 1-1/n (van Genuchten, 1980)}\cr
\code{01210}\tab{bimodel van Genuchten-Mualem model, with the contraint of m_i = 1-1/n_i (Durner, 1994)}\cr
\code{01310}\tab{trimodal van Genuchten-Mualem model, with the contraint of m_i = 1-1/n_i (Durner, 1994)}\cr
\code{02110}\tab{unimodel Kosugi 2 parametric-Mualem model (Kosugi, 1996)}\cr
\code{03110}\tab{unimodel Fredlund-Xing-Mualem model, with the contraint of m = 1-1/n (Fredlund D.G., and A. Xing, 1994)}\cr
}
}

\item{ivap.query}{
 \tabular{lll}{
\code{NULL}\tab{no isothermal vapour conductivity will be calculated with \code{Kvap}}\cr
\tab{Model type for isothermal vapour conductivity, see Details of function \code{KvapFun} for model codes}\cr
}
}
}
\details{
If the \code{shpmodel} code is supplemented by \code{FM}, e.g. \code{shpmodel = "01110"} as \code{shpmodel = "01110FM"}, the Framework-Model (Brunswick-Model) will be activated. Ensuring a water content of 0 at oven dryness (at a pressure head of h0 = 10^6.8 cm), continuous differentiability of the retention curve,
accounting for capillary and non-capillary water storage and conductivitey by introducing only one additional model parameter, and thepossibility 
to model a the often observed change in slope in the conductivity model in the medium pressure head range.

}
\value{ returns a \code{list} with calculations at specified \code{h}:
\item{theta}{calculated volumetric moisture content}
\item{Se}{calculated saturation}
\item{Scap}{effective saturation (of the capillary part if \code{FM} is specified}
\item{cap}{specific water capacity function (of the capillary part if \code{FM} is specified}
\item{psd}{pore size distribution (of the capillary part if \code{FM} is specified}
\item{Kh}{total hydraulic conductivity}
if \code{FM} specified, additionally:

\item{thetacap}{calculated volumetric moisture content of the capillary part}
\item{thetanc}{calculated volumetric moisture content of the non-capillary part}
\item{Snc}{effective saturation of the non-capillary part}
\item{Kcap}{hydraulic conductivity of the capillary}
\item{Knc}{hydraulic conductivity of the non-capillary}
\item{Kvap}{isothermal vapour conductivity}
\item{Krcap}{relative hydraulic conductivity of the capillary}
\item{Krnc}{relative hydraulic conductivity of the non-capillary}


}
\references{
\bold{Durner, W.}: Hydraulic conductivity estimation for soils with heterogeneous pore structurre, WRR, 30(2), 211-223, <doi:10.1029/93WR02676>, 1994.

\bold{Fredlund D.G., and Xing, A.}: Equations for the soil-water characteristic curve, Can. Geotech. J., 31:521-532, <doi:10.1139/t94-061>,1994.

\bold{Kosugi, K.}: Lognormal distribution model for unsaturated hydraulic properties, Water Reourc. Res., 32(9), 2697-2703, <doi:10.1029/96WR01776>, 1996.

\bold{Mualem, Y.}: New model for predicting hydraulic conductivity of unsaturated porous media. Water Resour. Res. 12(3): pp. 513-522. <doi:10.1029/WR012i003p00513>, 1976

\bold{Priesack, E. and Durner, W.}: Closed-Form Expression for the Multi-Modal Unsaturated Conductivity Function. Vadose Zone Journal 5:121-12, <doi:10.2136/vzj2005.0066>, 2006.

\bold{van Genuchten, M.T.}: Closed-form equation for predicting the hydraulic conductivity of unsaturated soils, Soil Sci Soc Am J, 44(5), 892-898, <doi:10.2136/sssaj1980.03615995004400050002x>, 1980.\cr

\bold{Weber, T.K.D., Durner, W., Streck, T., and Diamantopoulos, E.}:A modular framework for modelling unsaturated soil hydraulic properties over the full moisture range, WRR, in revision, 2019.
 }
\author{
Tobias KD Weber
}

\note{the function is used to assign a new function variable with a function which calculates the soil hydraulic properties 
according to sepcified \code{shpmodel} and model specified by \code{ivap.query}
}

\seealso{PC-Progress RETC
}
\examples{

# load measured data
data("shpdata1")
retdata <- shpdata1$LFH1$wrc[!is.na(shpdata1$LFH1$wrc[,1]),]
condata <- shpdata1$LFH1$hcc

# assign auxiliary variables
pF <- seq(-3, 6.8, length = 501)
h <- 10^pF

# assign list of parameters for the van Genuchten-Mualem model
parL <- list("p" = c("thr"= 0.05, "ths" = 0.45, "alf1" = 0.01, "n" = 2, "Ks" = 100, "tau" = .5),
 "psel" = c(1, 1, 0, 1, 1, 1),
 "plo" = c(0.001 , 0.2, 0.001, 1.1, 1, -2),
 "pup" = c(0.3, 0.95, 1, 10, 1e4, 10))

# calculate soil hydraulic property function values

shyp.L <- shypFun(parL$p, h, shpmodel = "01110", ivap.query = NULL)
wrc <- shyp.L$theta
hcc <- log10(shyp.L$Kh)

# # PLOT THE MEASURED WATER RETENTION CURVE
ticksatmin <- -2
tcllen <- 0.4
ticksat <- seq(ticksatmin,5,1)
pow <- ticksatmin:6

par(mfrow = c(1,2), tcl = tcllen)
plot(retdata, ylim = c(.0,.50), xlim = c(0, 6.8), ylab = "", xlab = "",
 col = "darkgrey", axes = FALSE, main = "Water Retention Curve", cex = 2)
lines(log10(h), wrc, col = "darkblue", lwd = 2)
legend("topright", c("observed", "simulated"),pch = c(1,NA),
 lty = c(NA,1), lwd = 2, bty = "n", cex = 1.3,col = c("darkgrey", "darkblue"))
axis(1, at = pow, labels = parse(text = paste('10^',(pow), sep = "")), tcl = tcllen)
axis(2, tcl = tcllen)
axis(4, labels = NA)
axis(3, labels = NA)
mtext("pressure head |h| [cm]", 1, cex = 1.2, line = 2.8)
mtext("vol. water content [ - ]", 2, cex = 1.2, line = 2.8)
box()

# PLOT THE MEASURED HYDRAULIC CONDUCTIVITY CURVE
plot(condata, ylim = c(-8,2), xlim = c(0, 6.8), ylab = "", xlab = "", col = "darkgrey",
 axes = FALSE, main = "Hydraulic Conductivity Curve", cex = 2)
lines(log10(h), hcc, col = "darkblue", lwd = 2)
legend("topright", c("observed", "simulated"), pch = c(1,NA),
 lty = c(NA,1), lwd = 2, bty = "n", cex = 1.3, col = c("darkgrey","darkblue"))
axis(1, at = pow, labels = parse(text = paste('10^',(pow), sep = "")), tcl = tcllen)
axis(2)
axis(4, labels = NA)
axis(3, labels = NA)
mtext("log10 K [cm/d]", 2, cex = 1.2, line = 2.8)
mtext("pressure head |h| [cm]",1 , cex = 1.2, line = 2.8)
box()
par(mfrow = c(1,1))


# HOW TO WRITE A MATER.IN FOR HYDRUS-1D

mater_out <- cbind(shyp.L[['theta']], h, shyp.L[['Kh']], abs(shyp.L[['cap']]))

materWriteFun <- function(mater_out.L, path = getwd(), sample) {

 # Function to write a Mater.in

 # ARGUMENTS

 # mater_outdata frame of 4 columns of calculated SHP values at h and length n. 
 # 1. Column: THETA
 # 2. Column: H(negative pressure heads)
 # 3. Column: K
 # 4. Column: C(positive)
 # path character specifying the path where the MATER.IN should be saved
 # sample optional chr for sample name: NULL = no name given

n <- dim(mater_out)[1]
sink(file.path(path, paste(sample, "MATER.IN", sep = "")))
cat(c("iCap", "\n", "1", "\n", "NTab", "\n", n, "\n"))
cat(c("\t","THETA", "\t\t\t","H","\t\t\t","K","\t\t\t","C"))
cat("\n")
 
write.table(format(mater_out, justify = "right"),
row.names = FALSE, col.names = FALSE, quote = FALSE)
sink()
 }
}
\keyword{ ~VGM}
