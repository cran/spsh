## Wrapper for all soil hydraulic property functions
shypFun <- function(p, h, shpmodel = "01110", ivap.query = NULL) {
          
          ### ASSIGN GENERAL 
          
          if(!is.na(grep("01110", shpmodel)[1]) == 1){
                    
                    # assign("shypFun", shypFun.01110)
                    shypFun.int <- shypFun.01110
          }
          
          
          if(!is.na(grep("01210", shpmodel)[1]) == 1){
                    
                    shypFun.int <- shypFun.01210
          }
          
          if(!is.na(grep("01310", shpmodel)[1]) == 1){
                    
                    shypFun.int <-  shypFun.01310
          }
          
          
          if(!is.na(grep("02110", shpmodel)[1]) == 1){
                    
                    shypFun.int <-  shypFun.02110
          }
          
          if(!is.na(grep("03110", shpmodel)[1]) == 1){
                    
                    shypFun.int <-  shypFun.03110
          }
          
          #### Simple calculation of soil hydraulic property functions without 
          #### the Framework model.
          
          if(is.na(grep("FM", shpmodel)[1] == 1)){
                    
                              shypL <- shypFun.int(p, h) 
                              kvap = KvapFun(p, por = p[2], retFun = NA, shypL$theta, model = ivap.query, pF = log10(h), output ="nolog10")
                              
                              return(list("theta" = shypL$theta,   
                                          "Scap" = shypL$Se,
                                          "cap" = shypL$cap, "psd" = shypL$psd,
                                          "Kh" = shypL$Kh + kvap, "Kcap" = shypL$Kh, "Kvap" = kvap
                              ))
                              
                              
                              
                              
                              
          }
          
          #### ADD SOMETHING IF FM IS ACTIVATED
          
          if(!is.na(grep("FM", shpmodel)[1] == 1)) {

                   {
                              h.meas <- abs(h)
                              h.support <- 10^seq(-3.11, 6.8, length = 302)
                              
                              h.calc <- sort(c(h.support, h.meas))
                              
                              index.h <- match(h.meas, h.calc)# %>%
                                        #unique(.) 
                              # Equations from Weber et al. (2018), submitted
                                       
                              thscnc <- p[1]            # [-]
                              thsc <- p[2]              # [-]
                              
                              # conductivity parameters
                              pcon <- p[(length(p)-4):length(p)]
                              
                              Ksc <- pcon[1]            # [cm d-1]
                              Ksnc <- pcon[3]           # [cm d-1]
                              afilm <- pcon[4]          # Set to 1.5, as suggested by Tokunaga et al. (2009), WRR
                              
                              h0 <- abs(10^pcon[5])     # Set to 6.8 (pF = log10(10^6.8)) )
                              #
                              ### SATURATION FUNCTION
                              #
                              #
                              
                              ## CAPILLARY PART, based on equation (13) 
                              #
                              u <- shypFun.int(p, h.calc)$Se
                              u0 = u[length(h.calc)]
                              
                              scap = (u-u0)/(1-u0);
                              
                              ## NON CAPILLARY part, based on equations (4-7)
                              #
                              snc <- sncFun(h.calc, scap)
                              
                              #
                              ### SPECIFIC WATER CAPACITY FUNCTION (Not in the manuscript)
                              #
                              #
                              
                              # Capillary part (NUMERICALLY)
                              #
                              capcap <- thsc * (diff(scap)/diff(h.calc))
                              nh <- length(capcap)
                              capcap <- c(capcap, capcap[nh-1]/2)                               # vector length adjustment
                              
                              # Non-capillary part (NUMERICALLY)
                              #
                              capnc <- thscnc * (diff(snc)/diff(h.calc))
                              nh <- length(capnc)
                              capnc <- c(capnc, capnc[nh-1]/2) 
                              
                              # Summation of both parts
                              cap <- capcap + capnc
                              
                              # pore size density after Durner 1994, WRR (valid only of the capillary part)
                              poredis = log(10) * h.calc * capcap
                              
                              #
                              ### SOIL WATER RETENTION CURVE, based on equation 1
                              #
                              #
                              thetacap = thsc * scap
                              thetanc = thscnc * snc
                              theta = thetacap + thetanc
                              
                              #
                              ### HYDRAULIC CONDUCTIVITY CURVE
                              #
                              #
                              
                              ## Capillary conductivity, based on equations (10-11)
                              #
                              Kcap <- shypFun.int(p, h.calc)$Kh
                              ## non-capillary conductivity similar to Peters (2013)
                              #
                              krnc = h0^( -afilm * (1-snc));
                              
                              # total non capillary conductivity, second summand in equation (9)
                              Knc <- Ksnc * krnc;
                              
                              ## isothermal vapour conductivity acc. to Sposito (2006), with varying permeability models
                              # calculation of isothermal vapor conductivity based on model MQ61 (refer to paper, Table A-1)
                              # equations A-1 to A-4 and Table A-1 model "MQ61
                              kvap = KvapFun(p, por = thscnc + thsc, retFun = NA, theta, model = ivap.query, pF=log10(h.calc), output ="nolog10")
                              
                              ## total unsat. conductivity based on equation (8)
                              Kh <-  Kcap + Knc + kvap
                              
                              ### Return list of all variables
                              return(list("theta" = theta[index.h], "thetacap" = thetacap[index.h], "thetanc" = thetanc[index.h],  
                                          "Scap" = scap[index.h], "Snc" = snc[index.h], 
                                          "capcap"=capcap[index.h], "capnc" = capnc[index.h],
                                          "cap" = cap[index.h], "psd" = poredis[index.h],
                                          "Kh" = Kh[index.h], "Kcap" = Kcap[index.h], "Knc" = Knc[index.h], "Kvap" = kvap[index.h],
                                          "krcap" = (Kcap/Ksc)[index.h], "krnc" = krnc[index.h]))
                    }
          }
}




