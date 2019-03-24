# returns the numerical solution to the Mualem 1976 integral
# Runs fast because of Riemann Integral
numMualem <- function(h, pcon, scap) {

          #
          ####  PURPOSE
          #
          #     Function to calculate the capillary bundle model, analytically
          #

          #
          #### ARGUMENTS 
          #
          #       h       vector of length l with pressure heads [cm]
          #       pcon    vector of conductivity model parameters, second entry has to be the parameter tau
          #       scap    vector of length l with non (!) rescaled saturation values for the capillary part
          
          #
          #### RETURNS vector of length l
          #
          #       mua_num     calculated capillary conductivity model 

          # general capacity conductivity parameters (Hoffmann-Riem et al., 1999)
          q <- 1
          r <- 2

          int1 = cumtrapz(rev(scap), rev(h^-q))
          int2 = max(int1)[1];
          # mua_num = rev((int1/int2)^r)
          
          # int1 = cumsum(rev(h^-q) * c(diff(rev(scap)[1:2])/2, diff(rev(scap))))
          # int2 = int1[length(int1)]
          mua_num = rev((int1/int2)^r)
          
          return(mua_num)
}


