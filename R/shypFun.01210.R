# RETC: constrained bimodal van Genuchten (1980), Durner (1994): constrained refers to substituting the shape parameter 'm' with '1-1/n'
# COND: Mualems Capillary Bundle Model (1976)
# Number of Parameters: 9

shypFun.01210 <- function(p, h){
          
          h <- abs(h)
          thr <- p[1]
          ths <- p[2]
          # first modality
          alf1 <- p[3]
          n1 <- p[4]
          w1 <- p[5]
          # second modality
          alf2 <- p[6]
          n2 <- p[7]
          # conductivity
          Ks <- p[8]
          tau <- p[9]
          
          m1 <- 1-1/n1
          m2 <- 1-1/n2
          w2 <- 1-w1
          
          Se1 = (1+(alf1*h)^n1)^-m1   
          Se2 = (1+(alf2*h)^n2)^-m2
          
          Se = w1*Se1 + w2*Se2  
          theta = thr + (ths-thr)*Se
          
          cap1 =      w1*alf1*m1*n1*(alf1*h)^(n1-1)*(1+(alf1*h)^(n1))^(-m1-1)
          cap2 =      w2*alf2*m2*n1*(alf2*h)^(n2-1)*(1+(alf2*h)^(n2))^(-m2-1)
          cap = (ths - thr) * (cap1+cap2)
          
          poredis = log(10)*h*cap
          
          Kh <- Ks * Se^tau *
                    ((w1 * alf1 * (1-(1-Se^(1/m1))^m1) +
                                w2 * alf2 * (1-(1-Se^(1/m2))^m2))/(w1*alf1+w2*alf2))^2 
          
          return(list("theta" = theta, "Se" = Se, "cap" = cap, "psd" = poredis, "Kh" = Kh))
}

