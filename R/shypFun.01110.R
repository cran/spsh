# RETC: constrained unimodal van Genuchten (1980): constrained refers to substituting the shape parameter 'm' with '1-1/n'
# COND: Mualems Capillary Bundle Model (1976)
# Number of Parameters: 6

shypFun.01110 <- function(p, h){
          
          h <- abs(h)
          thr <- p[1]
          ths <- p[2]
          alf1 <- p[3]
          n1 <- p[4]
          Ks <- p[5]
          tau <- p[6]
          
          m1 <- 1-1/n1
          Se <- (1+(alf1*h)^n1)^-m1        
          theta <- thr + (ths-thr) * Se  
          
          cap = -(ths-thr) * (n1*m1*alf1*(alf1*h)^(n1-1)*(1+(alf1*h)^(n1))^(-m1-1))
          
          poredis = log(10) * h * -cap
          Kh <-  Ks*Se^tau*(1-(1-Se^(1/m1))^m1)^2 
          
          return(list("theta" = theta, "Se" = Se, "cap" = cap, "psd" = poredis, "Kh" = Kh))
}

