sncFun <- function(h, scap){
          
          #
          ####  PURPOSE
          #
          #     Function to calculate the non-capillary saturation
          #
          
          #
          #### ARGUMENTS 
          #
          #       scap    vector of length l with non (!) rescaled saturation values for the capillary part
          
          #
          #### RETURNS vector of length l
          #
          #       snc     calculated non-capillary saturation values
          
          
          nh <- length(h)
          
          # solution to the integral
          # Snc_star_int =  pracma::cumtrapz(h,(1-scap)/h);
          # snc = 1-(Snc_star_int/Snc_star_int[nh]);
          # 
          result = tryCatch({
                Snc_star_int =  drop(pracma::cumtrapz(h,(1-scap)/h));
                1-(Snc_star_int/Snc_star_int[nh]);
          }, warning = function(w) {
                rep(1, nh)
          }, error = function(e) {
                rep(1, nh)
          }
          )
          # is.na(result) { result <- }
          return(result)
}
