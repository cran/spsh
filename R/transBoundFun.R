transBoundFun <- function(parL, shpmodel, weightmethod){
          
          
          #### AUXHILLARY VARIABLES 
          p    <- parL[[1]]           
          psel <- parL[[2]]
          plo  <- parL[[3]]
          pup  <- parL[[4]]
          
          pretransbase <- ptransbase <- rep(list(function(x)x), length(p))
          
          switch(gsub("FM", '',shpmodel),
                 "01110" = { 
                           # TRANFORMATION RULES OF PARAMETERS FOR van Genuchten-Mualem m = 1-1/n constrained MODEL, shpmodel = 01110
                           ptransfit    <- c(function(x)x, function(x)x,    log10,        function(x)log10(x-1),      log10     , function(x)x)
                           pretransfit  <- c(function(x)x, function(x)x, function(x)10^x, function(x)10^x+1    , function(x)10^x, function(x)x)

                 },
                 "01210" = {
                           ptransfit    <- c(function(x)x, function(x)x,    log10,        function(x)log10(x-1), function(x)x,    log10,        function(x)log10(x-1),      log10     , function(x)x)
                           pretransfit  <- c(function(x)x, function(x)x, function(x)10^x, function(x)10^x+1    , function(x)x, function(x)10^x, function(x)10^x+1    ,function(x)10^x, function(x)x)
                 },
                 "01310" = {  ptransfit    <- c(function(x)x, function(x)x,    log10,        function(x)log10(x-1),function(x)x,    log10,        function(x)log10(x-1), function(x)x,    log10,        function(x)log10(x-1),      log10     , function(x)x)
                              pretransfit  <- c(function(x)x, function(x)x, function(x)10^x, function(x)10^x+1    ,function(x)x, function(x)10^x, function(x)10^x+1    , function(x)x, function(x)10^x, function(x)10^x+1    ,function(x)10^x, function(x)x)
                 },
                 
                 "02110" = {
                           # TRANFORMATION RULES OF PARAMETERS FOR Kosugi-Mualem m = 1-1/n constrained MODEL, shpmodel = 01110
                           ptransfit    <- c(function(x)x, function(x)x,    log10,           log10           ,      log10     , function(x)x)
                           pretransfit  <- c(function(x)x, function(x)x, function(x)10^x, function(x)10^x    , function(x)10^x, function(x)x)
                           
                 },
                 "03110" = { 
                           # TRANFORMATION RULES OF PARAMETERS FOR Fredlund-Xing-Mualem m = 1-1/n constrained MODEL, shpmodel = 01110
                           ptransfit    <- c(function(x)x, function(x)x,    log10,        function(x)log10(x-1),      log10     , function(x)x)
                           pretransfit  <- c(function(x)x, function(x)x, function(x)10^x, function(x)10^x+1    , function(x)10^x, function(x)x)
                           
                 },
                 stop("Enter a valid code for a soil hydraulic property model")
          )
          
          
          #
          # query if FM selected and add extra parameters
          if(grepl("FM", shpmodel, fixed = TRUE)){
                    
                    ptransfit <- c(ptransfit, log10, function(x)x ,function(x)x)
                    pretransfit <- c(pretransfit,  function(x)10^x, function(x)x ,function(x)x)
          }


          if(!is.list(weightmethod)) {
                    
                    
                    if(weightmethod == "est1"){
                              
                              ptransfit   <-  c(ptransfit, log10, log10)
                              pretransfit <-  c(pretransfit, function(x)10^x, function(x)10^x)
                              
                              # check if user has provided additional parameters and bounds for the estimated standard deviation
                              # stopifnot(length(p) == ptransfit & length(psel) == ptransfit & length(plo) == ptransfit & length(pup) == ptransfit )
                              
                              # stopifnot(sum(sapply(parL, length))/length(parL) == length(ptransfit))
                              
                    }
          }
          
          ptransbase[which(psel == 1)]     <- ptransfit[which(psel == 1)] 
          pretransbase[which(psel == 1)]   <- pretransfit[which(psel == 1)] 
          
          # Modify the boundaries to exakt value for non selected parameters, only necessary if assignment in resFun exact to psel parameters
          plo[which(psel == 0)] <- unlist(p)[which(psel == 0)]
          pup[which(psel == 0)] <- unlist(p)[which(psel == 0)]
          
          return(list("parL" = list("p" = p, "psel" = psel, "plo" = plo, "pup" = pup), "transL" = list("ptrans" = ptransbase, "pretrans" = pretransbase )))
          
}
