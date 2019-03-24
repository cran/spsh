ptf.vG2FM <- function(x) {

      #
      ## Purpose parameter transfer function from van Genuchten Muaelem to Weber et al.
      #
      #
      
      #
      ## Arguments
      #
      #     x     num         vector with numerical parameter values
      #           1. theta r
      #           2. theta_s
      #           3. alf1
      #           4. n1
      #           5. Ks
      #           6. tau
      
      # Slope and Intercept
      M.t <- rbind("ptf.thr" = c(1.0024,0.03780),
      "ptf.ths" = c(0.9697, 0.00564),
      "log10Ks" = c(1.0280,-0.12968),
      "log10alf1" = c(0.8557,-0.29466),
      "log10n-1" = c(0.7178,0.18336),
      "tau" = c(0.8418, 0.61881))
      
      colnames(M.t) <- c("slope", "intercept")
      
      ptransfit    <- c(function(x)x, function(x)x,    log10,        function(x)log10(x-1),      log10     , function(x)x)
      pretransfit  <- c(function(x)x, function(x)x, function(x)10^x, function(x)10^x+1    , function(x)10^x, function(x)x)
      
      x_trans <- mapply(function(f, x) f(x), ptransfit, x)
      
      x.FM_trans <- c(M.t[,1]*x_trans + M.t[,2])
      
      x.FM <- mapply(function(f, x) f(x), pretransfit, x.FM_trans)
      
      w = .01
      
      x.FM.return <- c("thr" = x.FM[1],
                "ths" = (x.FM[2]-x.FM[1]),
                "alf1" = x.FM[3], 
                "n1" = x.FM[4],
                "Kcs" = x.FM[5]*(1-w),
                "tau" = x.FM[6],
                "Kncs" = x.FM[5]* w,
                "a" = 1.5,
                "h0" = 6.8)
      
return(x.FM.return)      
}
