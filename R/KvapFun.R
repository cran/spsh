KvapFun <- function(p, por = p[2], retFun = NA, theta = NA, model = "MQ61", Temp = 20, m = 3, pF = seq(-3, 7, length = 501), output = "log10",  ...) {
                    
                    # p, theta, h, por , Temp = 20, m = 3, model = "MQ61"
          
          if(is.null(model)){
                    return(0)
          }
          
          h = -10^pF
          # initial conditions
          Temp = 273.15 + Temp  # [ K ] temperature of experiment
          R = 8.314             # [ J mol^-1 kg^-1 ] 
          M = 0.018015          # [ kg mol^-1 ]
          roh_w = 1000*(1.0-7.37e-6*(Temp-273.15-4.0)^2.0+3.79E-8*(Temp-273.15-4.0)^3.0)   # [ kg m^-3]
          roh_vs = 0.001/Temp*exp(31.3716-6014.79/Temp-0.00792495*Temp);                   # [ kg m^-3]
          g = 9.81              # [ m s-^2 ]
          #tfac = 86400          # conversion factor [ m s^-1 ] --> [ m d^-1 ] 
          MGRT = M*g/(Temp*R)   # 
          Hr = exp(-abs(h*0.01)*MGRT)
          Da = 86400*2.14e-5 *(Temp/273.15)^2    # Diffusivity of water in air   (often denoted as Do or Da in literature)
          
          if(is.na(retFun) == FALSE) {
                    # calcualte air filled porosity 'eps(h)'
                    eps <- por - retFun(p, h)[['theta']]
                    # calcualte air filled porosity 'eps(h=100)'
                    eps100 <- por - retFun(p, h = 100)[['theta']]
          }
          
          else {eps = por - theta }
          
          switch(model,
                 B={          Do <- eps^2 * Da     
                 },
                 P={          Do <- 0.66 * eps * Da
                 },
                 MQ60={       Do <- (eps^2 * Da)/por^(2/3)
                 },
                 MQ61={       Do <- (eps^(10/3)*Da)/por^2
                 },
                 GS={         Do <- 1e-6*(eps*100)^3.36 * Da
                 },
                 L={          Do <- eps^(5/3) * Da
                 },
                 Xu={         Do <- eps^2.51/por^2 * Da
                 },
                 PMQ={        Do <- 0.66 * eps* (eps/por)^((12-m)/3)
                 },
                 # MC={       D <- (2*eps100^3 + 0.04 * eps100) * (eps/eps100)^(2+3/b)
                 # },         Find find a pedotransfer function for b
                 TPM={        X <- log10((2*eps100^3 + 0.04 * eps100)/por^2)/log10(eps100/por)  
                              Do <- por^2 * (eps/por)^X
                  },
                 TPEM={       X <- 2 + log10(eps100^0.25)/log10(eps100/por)
                              Do <- por^2 * (eps/por) ^ X
                  },
                 stop("Enter something that switches the model in KvapFun")
          )
          
          # Calculate the isothermal vapour conductivity in cm/d
          Kvap = 100 * roh_vs/roh_w * Do * MGRT * Hr
          
          # return
          if(output == "log10") {
                    return(log10(Kvap)) 
          }
          else {return(Kvap)}
}


# Kvap <- KvapFun(unlist(p),por = p[2], theta = output$theta)

