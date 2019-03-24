# RETC: constrained trimodal van Genuchten (1980), (Durner 1994) (c.f. Priesack and Durner for analytical solution): constrained refers to substituting the shape parameter 'm' with '1-1/n'
# COND: Mualems Capillary Bundle Model (1976)
# Number of Parameters: 12

shypFun.01310 <- function(p, h){
          
          h <- abs(h)
          thr <- p[1]
          ths <- p[2]
          
          alf1 <- p[3]
          n1 <- p[4]
          w1 <- p[5]
          
          alf2 <- p[6]
          n2 <- p[7]
          w2 <- p[8]
          
          alf3 <- p[9]
          n3 <- p[10]
          
          Ks <- p[11]
          tau <- p[12]
          q <- 1
          r <- 2

          # computed parameters
          ah1 <- alf1*h
          ah2 <- alf2*h
          ah3 <- alf3*h
          
          m1 <- 1-1/n1
          m2 <- 1-1/n2
          m3 <- 1-1/n3
          
          w3 <- 1-w1-w2
          
          # capillary storage
          u1 = (1+ah1^n1)^-m1   
          u2 = (1+ah2^n2)^-m2
          u3 = (1+ah3^n3)^-m3
          scap = w1*u1 + w2*u2  + w3*u3
          # water content
          theta = thr + (ths-thr)*scap
          
          # capacity function
          cap1 =      w1*alf1*m1*n1*ah1^(n1-1)*(1+ah1^n1)^(-m1-1)
          cap2 =      w2*alf2*m2*n2*ah2^(n2-1)*(1+ah2^n2)^(-m2-1)
          cap3 =      w3*alf3*m3*n3*ah3^(n3-1)*(1+ah3^n3)^(-m3-1)
          cap = (ths - thr) * (cap1+cap2+cap3)
          
          # pore distribution
          poredis = log(10)*h*cap
          
          # conductivity
          fs <- w1 * alf1 * (1-u1^(1/m1))^m1 + 
                    w2 * alf2 * (1-u2^(1/m2))^m2 + 
                    w3 * alf3 * (1-u3^(1/m3))^m3
          
          f0 <- w1*alf1 + w2*alf2 + w3*alf3
          
          mua <- 1-fs/f0
          Kh <- Ks * scap^tau * mua^r

          return(list("theta"=theta, "Se" = scap, "cap" = cap, "psd" = poredis, "Kh" = Kh))
}