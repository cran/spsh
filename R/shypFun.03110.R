# RETC: FFredlund D.G., and A. Xing. 1994. Equations for the soil-water characteristic curve. Can. Geotech. J. 31:521-532.
# COND: Mualems Capillary Bundle Model (1976)
# Number of Parameters: 6
# No conductivity calculation. This can be done through the shypFun.

shypFun.03110 <- function(p, h){
          
          h <- abs(h)
          thr <- p[1];
          ths <- p[2];
          alf1 <- p[3];           # cm^-1
          n1 <- p[4];
          Ks <- p[5]
          tau <- p[6]
          
          m1 <- 1-1/n1       # constrained case

          Se <- (log(exp(1)+(alf1*h)^n1))^-m1
          theta <- thr + (ths-thr) * Se  
          
          cap_num = diff(theta)/diff(h)
          cap <- c(cap_num[1],cap_num)
          
          poredis <- log(10) * h * -cap

          Kh <- Ks * Se^tau * numMualem(h, pcon = tau, scap = Se)
          
          return(list("theta" = theta, "Se" = Se, "cap" = cap, "psd" = poredis,
                      "Kh" = Kh))
}
