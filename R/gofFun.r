gofFun <- function(phat, shpmodel = "01110", retdata = NULL, condata = NULL,
                                        weight, psel, ivap.query = NULL, hclip.query = FALSE){
          
          # initialise two lists
          th <- vector("list")
          logKh <- vector("list")
          combined <- vector("list")                    # multi-objective goodness-of-fit
          
          #  - check for trimodal soil hydraulic property model. Included for future functionality
          ifelse(shpmodel%in%c("01310"), trim.query <- TRUE, trim.query <- FALSE)
          
          # calculate vectors of residuals
          res.L <- resFun(phat, shpmodel, retdata, condata, pretrans = NULL, weight, method = "res", trim.query, ivap.query, hclip.query) # residuals of THETA and log10Kh
          
          n <- sum(psel)                             # number of estimated SHP model parameters
          
          # 1. calculate for th
          if(!is.null(condata)){
                    th$res <- res.L[1:dim(retdata)[1]]       # Residuals of THETA
                    r <- th$res
                    m = length(r);
                    th$me   = mean(r);                       # mean (weighted) error
                    th$mae  = mean(abs(r));                  # mean average (weighted) error
                    th$mse  = sum(r^2)/m                     # mean squared (weighted) error
                    th$rss  = t(r)%*%r;                      # sum of squared (weighted) residuals
                    th$rmse = sqrt(th$rss/m);                # root mean squared (weighted) error
                    th$logLik = log(th$rss/m)
                    th$AIC  = -2 * th$logLik + 2 * n;
                    th$AICc = ifelse((m-n-1) <= 0,NaN, th$AIC + 2*n*(n+1)/(m-n-1));
                    th$BIC  = -2 * th$logLik+n*log(m);
                    th$m    = m;                             # Number of THETA-H-pairs
          }
          
          # 2. calculate for Kh
          if(!is.null(condata)){
                    logKh$res <- res.L[(dim(retdata)[1]+1) : (dim(retdata)[1]+dim(condata)[1]) ]             # Residuals of log10Kh
                    r = logKh$res
                    m = length(r);
                    logKh$me   = mean(r);                    # mean (weighted) error
                    logKh$mae  = mean(abs(r));               # mean average (weighted) error
                    logKh$mse  = sum(r^2)/m;                 # mean squared (weighted) error
                    logKh$rss  = t(r)%*%r;                   # sum of squared (weighted) residuals
                    logKh$rmse = sqrt(logKh$rss/m);          # root mean squared (weighted) error
                    logKh$logLik = log(logKh$rss/m)
                    logKh$AIC  = -2*logKh$logLik + 2*n;
                    logKh$AICc = ifelse((m-n-1) <= 0, NaN, logKh$AIC + 2*n*(n+1)/(m-n-1));
                    logKh$BIC  = -2 * logKh$logLik + n*log(m);
                    logKh$m    = m;                          # Number of LOG10K-H-pairs
          }
          
          
          # returns
          if(!is.null(retdata) & !is.null(condata)){
                    combined$AIC  <- -2 * (log(th$rss/m) + log(logKh$rss/m))  + 2 * n;
                    combined$AICc <- ifelse((m-n-1)<=0,NaN,combined$AIC + 2 * n * (n+1) / (m-n-1));
                    combined$BIC  <- -2 * (log(th$rss/m) + log(logKh$rss/m)) + n * log(m)
                    
                    return(list("th" = th, "log10Kh" = logKh, "combined" = combined))
          }
          
          if(!is.null(retdata) & is.null(condata)){
                    return(list("th" = th, "log10Kh" = NULL))
          }
          
          if(is.null(retdata) & !is.null(condata)){
                    return(list("th" = NULL, "log10Kh" = logKh))
          }
}



