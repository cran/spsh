# RETC: 2 parametric Kosugi function (Kosugi, 1996)
# COND: Mualems Capillary Bundle Model (1976)
# Number of Parameters: 6

shypFun.02110 <- function(p, h){
          
          h <- abs(h)
          thr <- p[1];
          ths <- p[2];
          hm <- p[3];           # cm^-1
          sigma <- p[4];
          Ks <- p[5]
          tau <- p[6]

          Se <- 0.5 * erfc(log(h/hm)/(sqrt(2)*sigma))
          
          theta <- thr + (ths-thr) * Se  
          
          cap <- (ths-thr)/(sqrt(2*pi)*sigma*h)*exp(log(h/hm)^2/(2*sigma^2))
          
          poredis <- log(10) * h * -cap
          Kh <-   Ks * Se^tau * (0.5 * erfc(log(h/hm)/(sqrt(2)*sigma) + sigma/sqrt(2) )  )^2
          
          return(list("theta" = theta, "Se" = Se, "cap" = cap, "psd" = poredis, "Kh" = Kh))
}
