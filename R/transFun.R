transFun <- function(par.vec, trans.L) {
          
          #
          ####  PURPOSE
          #
          #     Function to transform and back-transform the parameters (of e.g.  soil hydraulic property functions),
          #         this enables more robust estimation by reducing the correlation structure of the parameters
          
          #
          #### ARGUMENTS 
          #
          #       par.vec     num       vector of length l of parameters to be transformed
          #       trans.L     function  list of transformation functions of same length par.vec
          
          #
          #### RETURNS vector of length l
          #
          #       p.transform num       vector of length l with transformed parameters.
          
          p.transform <- mapply(function(f, x) f(x), trans.L, par.vec)
          
          return(p.transform)
          
}