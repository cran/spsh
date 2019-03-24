weightFun <- function(weightmethod = "fix1", retdata, condata, parL = NA){
          
      if(weightmethod == "est1"){
            {stopifnot(!is.na(parL))}
      }
      
          if(weightmethod == "range"){
                    weight <- list("wth" = max(retdata[,2]) - min(retdata[,2]), "wKh" = max(condata[,2]) - min(condata[,2]))
          }
          
          # NO WEIGHTING
          if(weightmethod == "none"){
                    weight <- list("wth" = 1, "wKh" = 1)
          }
          
          # MEASUREMENTS PREDEFINED DUE TO A FIXED MEASUREMENT UNCERTAINTY Peters and Durner (2008)
          if(weightmethod == "fix1"){
                    weight <- list("wth" = 0.0025, "wKh" = 1)
          }
          
          # MEASUREMENTS PREDEFINED DUE TO A FIXED MEASUREMENT UNCERTAINTY
          if(weightmethod == "fix2"){
                    weight <- list("wth" = c(rep(0.05, dim(retdata)[1])), "wKh" = c(rep(1,dim(condata)[1])))
          }
          
          # PARAMETERS ESTIMATED
          if(weightmethod == "est1"){
                    
                    p <- parL[[1]]
                    psel <- parL[[2]]
                    plo <- parL[[3]]
                    pup <- parL[[4]]
                    
                    weight <- NULL
                    # add 
                    pstat    <- c(sdth = 0.1 , sdKh = 0.1)
                    pstatsel <- c(   1       ,    1   )
                    pstatlo  <- c(   1e-3    ,    1e-3   )
                    pstatup  <- c(   1e0     ,    1e1    )
                    
                    p    <- c(p  , pstat   )
                    psel <- c(psel, pstatsel)
                    
                    plo <- c(plo, pstatlo)
                    pup <- c(pup, pstatup) 
                    
                    plo[which(psel == 0)] <- p[which(psel == 0)]
                    pup[which(psel == 0)] <- p[which(psel == 0)]
                    
                    # ptrans   <- c(ptrans, function(x)x, function(x)x)
                    # pretrans <- ptrans  
                    # 
                    # ptransfit   <- c(ptrans, log10, log10)
                    # pretransfit <- c(pretrans, function(x)10^x, function(x)10^x)
                    # 
                    # ptrans[which(psel == 1)]   <- ptransfit[which(psel == 1)] 
                    # pretrans[which(psel == 1)] <- pretransfit[which(psel == 1)] 
                    
          }

          if(weightmethod != "est1"){ 
                    return(weight)   
          } else {
                    return(list("parL" = list("p" = p, "psel" = psel, "plo" = plo, "pup" = pup)))
                    }
}


