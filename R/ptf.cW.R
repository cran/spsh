ptf.cW <- function(CLAY, SAND, BD, OC) {
          
          #
          ####  PURPOSE
          #
          #     Function to predict soil hydraulic property model parameters for the constrained van Genuchten Mualem function
          #         function based on Weynant et al. (2009) and corrections by Weiherm?ller
          
          #
          #### ARGUMENTS 
          #
          #         CLAY      num       clay fraction [0, 100] [% weight]
          #         SAND      num       sand fraction [0, 100] [% weight]
          #         BD        num       bulk density of the soil [g cm-3] 
          #         OC        num       organic carbon content [0, 100][% weight]
          
          #
          #### RETURNS list of length l
          #
          #       cW          num       list of 6 model parameters, each list of same length as the arguments
          
          
          # Internal tuning parameters
          par <- c(0.0, 0.0, 0.0, 0.6355365504, 0.0012807613, -0.1630910071,
                   -4.3003420026, -0.0096851210, 0.0137763939, 0.0, -0.0992391618,
                   -1.0846069190, -0.0235654376, -0.0085043044, 0.0001369944, 1.9581584626,
                   0.0, 0.0308449531, -0.6141610130, -0.1565557100, -1.8641969563,
                   -0.1317363417, 0.0067082506, 0.0, 0.0)
          
          cW <- list()
          
          cW$thr <- par[1] + par[2] * CLAY + par[3] * OC
          cW$ths <- par[4] + par[5] * CLAY + par[6] * BD
          lalf1 <- par[7] + par[8] * CLAY + par[9] * SAND + par[10] * BD + par[11] * OC
          ln1 <- par[12] + par[13] * CLAY + par[14] * SAND + par[15] * SAND^2
          
          cW$alf1 <- exp(lalf1)
          cW$n1 <- exp(ln1)+1
          
          lK <- par[16] + par[17] * CLAY + par[18] * SAND + par[19] * BD + 
                    par[20] * OC
          cW$K0 <- exp(lK)
          cW$tau = par[21] + par[22] * CLAY + par[23] * SAND + par[24] * BD + par[25] * OC
          
          return(cW)
}
