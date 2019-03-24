# 3 # Residual Functions  -----------------------------------------------------
# normal distributed residuals without weights: 
resFun <- function(p, shpmodel = "01110", retdata = NULL, condata = NULL, pretrans = NULL, 
                   weight = NULL, method = "rss", 
                   trim.query = FALSE, ivap.query = NULL, hclip.query = FALSE, parL = NULL) {
          #
          ####      ARGUMENTS
          #
          #         p         - vector of num transformed soil hydraulic model parameters
          #         retdata   - data.frame of num observations of retention data first column: pressure heads, second column: volumetric water content [-]
          #         condata   - data.frame of num observations of hydraulic conductivity data, data.frame, first column: pressure heads, second column: log10K [cm/d]
          #         pretrans  - list of functions with same length as parameter vector to back-transform the parameters.
          #         weight    - List of 2 vectors with weights (inverse of variance, describing e.g. the measurement erro)  
          #         method    - chr character specifying the adopted fitting method.
          #         trim      - trimodal query
          
          # assigning fixed paraemeters to the parameter vector
          
          if(!is.null(pretrans)){ 
                    p <- transFun(p, pretrans)
          }
          
          if(!is.null(parL)){ 
                    parL$p[which(parL$psel == 1)] <- p
                    p <- parL$p
          }
          
          # assigning the weights if they are estimated scalars and remove the parameters from the vector
          if(is.null(weight)){
                    np <- length(p)
                    sdth <- p[np-1];
                    sdKh <- p[np];
                    p <- p[-c(np-1, np)]
          }
          
          
          # query if trimodal function is used: if TRUE, the sum of weights is calculated, to verify sum == 1
          if(trim.query == TRUE){
                    w1 <- p[5]
                    w2 <- p[8]
                    if((1-w1-w2) <= 0) {
                              
                              if(method == "res") {
                                        return(rep(1000, dim(retdata)[1] + dim(condata)[1]))
                              } 
                              
                              if(method == "rss") {
                                        return(1e32)
                              }
                              
                              if(method == "-2logLik")
                                        return(1e32)
                    }
          }      
          
          # calculation of weights if they are fix
          if(!is.null(weight) & method == "rss" ) {
                    sdth <- 1/weight[[1]];
                    sdKh <- 1/weight[[2]];
          }
          
          if(!is.null(weight) & method == "res" ) {
                    sdth <- 1/weight[[1]];
                    sdKh <- 1/weight[[2]];
          }
          
          if(!is.null(weight) & method == "-2logLik" ) {
                    sdth <- weight[[1]];
                    sdKh <- weight[[2]];
          }
          
          # CALCULATE RESIDUALS
          yth <- retdata[,2]
          yhatth <- shypFun(p, h = 10^retdata[,1], shpmodel, ivap.query)[['theta']] 
          
          yKh <- condata[,2] 
          yhatKh <- log10(shypFun(p, h = 10^condata[,1], shpmodel, ivap.query)[['Kh']])
          
          # if(!isFALSE(ivap.query)) {
          #           if(!isTRUE(ivap.query)) {ivap.query <- "MQ61"}
          #           ifelse(ivap.query%in%c("MC", "TPM","TPEM"), retFun.query <- shypFun, retFun.query <- NA)
          #           yhatKh <- log10(10^yhatKh + KvapFun(p, por = p[2], retFun = retFun.query, theta = yhatth, model = ivap.query, Temp = 20, m = 3, pF = 10^condata[,1], output = "nolog10"))
          # }
          
          # OUTPUT
          if(method == "rss") {
                    res_th <- sdth * (yth - yhatth)
                    res_Kh <- sdKh * (yKh - yhatKh)
                    
                    rss <- drop(crossprod(res_th) + crossprod(res_Kh))
                    if(is.finite(rss)==TRUE) return(rss)
                    else(return(10^15))
          }
          
          if(method == "-2logLik"){
                    
                    # normal distribution
                    loglikth <- logLikFun.norm(yth, yhatth, sdth)
                    loglikKh <- logLikFun.norm(yKh, yhatKh, sdKh)
                    
                    if(sum(is.finite(c(loglikth, loglikKh))) == 2){
                              return(-2 * (loglikth + loglikKh))
                    }
                    if(sum(is.finite(c(loglikth, loglikKh))) != 2){
                              return(10^10)
                    }
          }
          
          if(method == "res") {
                    res_th <- sdth * (yth - yhatth)
                    res_Kh <- sdKh * (yKh - yhatKh)
                    res <- c(res_th,res_Kh)
                    if(sum(!is.finite(res)) == 0) return(res)
                    else(return(rep(10^6, length(res))))
          }
} 
