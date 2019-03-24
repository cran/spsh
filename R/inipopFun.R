# Function to generate an inital population on a latin hypercube


inipopFun <- function(p, psel, plo, pup, trans.L, Npop = NA) {
          
          # transform
          plotrans <- transFun(plo, trans.L)
          puptrans <- transFun(pup, trans.L)   
          
          # prepare "upscaling" from uniform distribution [0,1] to transformed parameter bounds [lower,upper]
          ran   <- puptrans - plotrans
          
          # Set population for lhs
          Npop  <- ifelse(is.na(Npop), sum(psel)*10, ceiling(Npop))
          
          # create population
          inipop0 <- optimumLHS(n = Npop, k = length(p), maxSweeps = 20, eps = .1, verbose = FALSE)
          
          # upscale to parameters | This works, since plo was set to the fixed parameter in p
          inipop  <- t(apply(inipop0 , 1, function(x) x*ran + plotrans))
          
          # return
          return(inipop)
}
