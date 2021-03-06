#' Measured soil hydraulic property data
#'
#' Data from evaporation experiments using the UMS HYPROP device
#' on two soils with different textures 
#' Data ist reported in [cm3 cm-3] 
#'
#' @docType data
#'
#' @usage data(shpdata1)
#'
#' @format An object of class /code{"list"}.
#'
#' @keywords datasets
#'
#' @references Kettcheson, Price, Weber (2018):  Initial variability, evolution and model parameterization 
#' of the soil hydrophysical properties of a reclaimed oil sands watershed and constructed wetland, in revision.
#'
#' @source none
#'
#' @examples
#' 
#' data(shpdata1)
#' ticksatmin <- -2
#' tcllen <- 0.4
#' ticksat   <- seq(ticksatmin,5,1) 
#' pow <- ticksatmin:5    
#'
#' TS.wrc <- shpdata1$TS1$wrc; TS.hcc <- shpdata1$LFH1$wrc
#' TS.hcc <- shpdata1$TS1$hcc; LFH.hcc <- shpdata1$LFH1$hcc
#'                                                                                            
#' # PLOT THE MEASURED WATER RETENTION CURVE
#' plot(TS.wrc[,1:2]  , ylim = c(.1,.50), xlim = c(0,4), ylab = "", xlab = "", 
#' col = "darkgrey", axes=FALSE, main = "Measured Water Retention Curve")
#' points(TS.wrc[,1:2], pch = 4,col = "darkblue")
#' legend("topright", c("TS1", "LFH1"), pch = c(1,4), bty = "n", cex=1.2, 
#' col = c("darkgrey","darkblue"))
#' axis(1, at = pow, labels = parse(text = paste('10^',(pow), sep = "")),tcl = tcllen)  
#' axis(2, tcl = tcllen)
#' axis(4, labels = NA, tcl = tcllen)
#' axis(3, labels = NA, tcl = tcllen)  
#' mtext("pressure head |h| [cm]", 1, cex = 1.2, line = 2.8)
#' mtext("vol. water content [ - ]", 2, cex = 1.2, line = 2.8)
#' box()
#'
#' # PLOT THE MEASURED HYDRAULIC CONDUCTIVITY CURVE
#'
#' plot(TS.hcc, ylim = c(-8,2), xlim = c(0,4), ylab = "", xlab = "", col = "darkgrey", 
#' axes = FALSE, main = "Measured Hydraulic Conductivity Curve")
#' points(TS.hcc, pch = 4, col = "darkblue")
#'
#' legend("topright", c("TS1", "LFH1"), pch = c(1,4), bty = "n", cex = 1.2, 
#' col = c("darkgrey","darkblue"))
#'
#' axis(1, at = pow, labels = parse(text = paste('10^',(pow), sep = "")), tcl = tcllen)  
#' axis(2, tcl = tcllen)
#' axis(4, labels = NA, tcl = tcllen)
#' axis(3, labels = NA, tcl = tcllen) 
#'
#' mtext("log10 K [cm/d]", 2, cex = 1.2, line = 2.8)
#' mtext("pressure head |h| [cm]",1 , cex = 1.2, line = 2.8)
#' box()
"shpdata1"
