
logLikFun.norm <- function(y, yhat, sigma){
          
          #
          #### ARGUMENTS
          #
          #   y       num     vector of observed quantity
          #   yhat    num     vector of model predicuted quantity
          #   sigma   num     standard deviation of the model residuals (y-yhat)
          #
          
          #
          #### PURPOSE
          #
          #   calculate the log-Likelihood value 
          #
          #### ASSUMPTIONS
          #
          #   1)  Residuals follow a normal distribution and are standardized to N~(0,1^2)
          #   2)  residuals are independent
          #
          #   Function can account heteroscedastic and homoscedastic residuals
          #   sigma   num     standard deviation of the model residuals (y-yhat)
          #
          
          
          #
          #### RETURNS
          #
          #   loglik_norm   num   skalar log-Likelihood value 
          
          N <- length(y)
          
          eta <- (y-yhat)/sigma
          
          loglik_norm_sum <- -2*(-N/2*log(2*pi) -N*log(sigma)- 1/2 * sum(eta^2))
          
          return(loglik_norm_sum)
          
}
