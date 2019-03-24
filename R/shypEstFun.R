shypEstFun <- function(shpmodel = "01110", parL, retdata, condata, 
                       ivap = NULL, hclip = FALSE, 
                       weightmethod = "none", LikModel = "rss",  ALG = "DE", 
                       set.itermax = 200, ALGoptions = NULL, lhs.query = FALSE){
          
          # shypEstFun <- function(likmod ="rss", shpmodel = "01110", retdata, condata, parL, weightmethod, ALG = "DE", set.itermax = 200){
          

          #
          ####      ARGUMENTS
          #
          #         shpmodel            chr       Soil hydrailic property model, defaults to "01110" van Genuchten-Mualem constrained 1-1/m
          #
          #         parL                list      of 4 vectors with same length
          #                                                 1) p    num        vector of initial parameter estimates | the value for the fixed parameters
          #                                                 2) psel num        vector with 0s and 1s specifiying estimated parameters | 0 = no estimation, 1 = estimation
          #                                                 3) plo  num        vector of lower parameter boundaries
          #                                                 4) pup  num        vector of upper parameter boundaries
          #
          #         retdata             data.frame       observations of retention data, data.frame, 
          #                                                 first column     num   pressure heads in pF scale [cm], 
          #                                                 second column    num   volumetric water content [-]
          #
          #         condata             data.frame       observations of hydraulic conductivity data, 
          #                                                 first column     num   pressure heads [cm], 
          #                                                 second column    num   conductivity values, log10 transfromed | log10K [cm/d]
          #
          #         ivap                logical   FALSE (default) isothermal vapour conductivity (ivap) is not considered.  
          #                                       TRUE ivap is considered, defaults to Millington-Quirk (1961)
          #                             chr       isothermal vapour conductivity as described in Saito et al. (2006) is considered with a given specified model 
          #
          #         hclip               logical   Select of hclip correction method for the conductivity data after Iden et al. (2015) will be done, defaults to FALSE 
          #                                       (for future functionality, only)
          #
          #         weightmethod        chr       Selection of weighting method for the inverse modelling, defaults to none
          #                                                 range     chr       the data range of data values is considered (cf. Peters ad Durner, 2008)
          #                                                 fix1      chr       
          #                                                 fit       chr       the variance of the selected data classes is estimated, assuming homoscedastic variance.
          #                                                 none      chr       no weighting is considered.
          #                                                 est 1     chr       the standard deviation is estimated as a nuisance parameter
          #                             list
          #
          #         LikModel            chr       likelihood model, defaults to "rss"
          #                                                 rss       chr       minimisation of the  sum of squared (weighted) residuals
          #                                                 loglik    chr       minimisation of the log-likelihood
          #
          #                                                 res       chr       returns vector of (weighted) residuals                                  
          #         ALG              chr       Specification of ptimisation/sampling algorithm (for future functionality).
          #                                                 DE        chr       Differential Evolution | DEoptim package
          #
          #         set.itermax         num       number of maximum iterations            list      later: fit.control control parameters; see DEoptim.control.          
          #         
          #         lhs.query           chr       FALSE or TRUE: asks if a Latin Hypercube sample should be drawn for initial population for DEoptim
          
          #
          ####      INPUT CHECKS
          #                   Stop function execution if the following checks are negative
          
          
          #  - Is shpmodel ill-specified?
          stopifnot(shpmodel%in%c("01110","01110FM",
                                  "01210",
                                  "01310",  
                                  "02110", "02110FM",
                                  "03110", "03110FM"))   
          
          
          #  - Is input data not of class data.frame?
          stopifnot(is.data.frame(retdata), is.data.frame(retdata))
          
          #  - Is likmod ill-specified?
          stopifnot(LikModel%in%c("rss", "res", "-2logLik"))
          
          #  - Is likmod ill-specified?
          stopifnot(ALG%in%c("DE", "modMCMC"))
          
          
          #  - Is weightmethod ill-specified?
          stopifnot(weightmethod%in%c("range", "fix1", "fix2", "est1", "none")| is.list(weightmethod) )
          
          #  - Is parL ill-specified? 
          stopifnot(is.list(parL))
          
          #  - Does parL contain vectors of different lengths?
          stopifnot(sum(sapply(parL, length))/length(parL) == length(parL[[1]]))
          
          #  - parL is not a list of four 
          stopifnot(length(parL) == 4)
          
          #  - check for trimodal soil hydraulic property model. Included for future functionality
          ifelse(shpmodel%in%c("01310"), trim.query <- TRUE, trim.query <- FALSE)
          
          #  - check if ALGoptions exists and if not create it
          if(is.null(ALGoptions)){ALGoptions <- list()}
          
          #  - assigning some options
          ALGoptions$itermax <- set.itermax
          if(!"storepopfrom" %in% names(ALGoptions)){ALGoptions$storepopfrom <- ALGoptions$itermax+1}
          
          #
          ####       ASSIGN CORRECT FUNCTION FOR resFun to operate properly
          #
          
          #
          ####       INVERSE MODELLING ----
          #
          
          
          #
          ###        WEIGHTING METHOD   
          #
          
          if(is.list(weightmethod)){            # USER DEFINED WEIGHTS
                    
                    weight = weightmethod
                    
                    # Ensure that weight is either a skalar, or given as the same length as number of observations in
                    #  retdata
                    stopifnot(length(weight[[1]]) == 1 | length(weight[[1]]) == dim(retdata)[1]) 
                    #  condata
                    stopifnot(length(weight[[2]]) == 1 | length(weight[[2]]) == dim(condata)[1]) 
                    
          } else {
                    
                    if(weightmethod != "est1"){ # FIXED SKALAR WEIGHTS
                              weight <- weightFun(weightmethod, retdata, condata)
                              
                    } else {  
                              weight <- NULL    # FITTED SKALAR SD of THETA and SD of log10K
                              parL <- weightFun(weightmethod, retdata, condata, parL = parL)$parL
                    }  
          }
          
          #
          #####      TRANSFORM BOUNDARIES AND GENERATE INITIAL POPULATION, HYDRAULIC PROPERTY MODEL SPECIFIC
          #
          
          # get transform parameters.
          returnL <- transBoundFun(parL, shpmodel, weightmethod)
          parL <- returnL$parL
          transL <- returnL$transL
          
          # get transform parameters.
          p <- parL$p 
          plo_bd <- parL$plo 
          pup_bd <- parL$pup 
          psel <- parL$psel
          
          
          # initial population on transformed parameters.
          if(ALG == "DE"){
                    if(lhs.query == TRUE) {
                              ALGoptions$initialpop <- inipopFun(p = p, psel = psel, plo = plo_bd , pup = pup_bd, trans.L = transL$ptrans, Npop = NA)          
                    }else{
                              if(!"initialpop"%in%names(ALGoptions)) ALGoptions$initialpop <- NULL}
          }

          #
          ###       FITTING 
          #
          
          switch(ALG,
                 DE = {out <- DEoptim(resFun, 
                  # model specific
                  lower = transFun(plo_bd, transL$ptrans),
                  upper = transFun(pup_bd, transL$ptrans),
                  shpmodel = shpmodel,
                  retdata = retdata, 
                  condata = condata, 
                  pretrans = returnL$transL$pretrans, 
                  method = LikModel,
                  weight = weight,
                  trim.query = trim.query,
                  ivap.query = ivap,
                  
                  # algorithm specific
                  # DEoptim.control(VTR = -Inf, 
                  #                 strategy = 2, 
                  #                 bs = FALSE, 
                  #                 NP = sum(parL$psel)*10,
                  #                 itermax = set.itermax, 
                  #                 CR = .8, 
                  #                 F = 0.8, 
                  #                 trace = set.itermax/10, 
                  #                 initialpop = inipop.m,
                  #                 storepopfrom = set.itermax + 1,
                  #                 storepopfreq = 200, 
                  #                 p = 0.2, 
                  #                 c = 0, 
                  #                 reltol = 1e-8,
                  #                 parallelType = 0,  # 1
                  #                 cluster = NULL, 
                  #                 packages = c(), 
                  #                 parVar = c(), 
                  #                 foreachArgs = list())
                  DEoptim.control(
                            VTR  = ifelse("VTR" %in% names(ALGoptions), ALGoptions$VTR, -Inf)
                            , strategy = ifelse("strategy" %in% names(ALGoptions), ALGoptions$strategy, 2)
                            , bs = ifelse("bs" %in% names(ALGoptions), ALGoptions$bs, FALSE)
                            , NP = ifelse("NP" %in% names(ALGoptions), ALGoptions$NP, NA)
                            , itermax = ALGoptions$itermax
                            , CR = ifelse("CR" %in% names(ALGoptions), ALGoptions$CR, .5)
                            , F = ifelse("F" %in% names(ALGoptions), ALGoptions$F, .8)
                            , trace = ifelse("trace" %in% names(ALGoptions), ALGoptions$trace, TRUE)
                            , initialpop = if("initialpop" %in% names(ALGoptions)){ALGoptions$initialpop}
                            , storepopfrom = ALGoptions$storepopfrom
                            , storepopfreq = ifelse("storepopfreq" %in% names(ALGoptions), ALGoptions$storepopfreq, 1)
                            , p = ifelse("p" %in% names(ALGoptions), ALGoptions$p, .2)
                            , c = ifelse("c" %in% names(ALGoptions), ALGoptions$c, 0)
                            , reltol = ifelse("reltol" %in% names(ALGoptions), ALGoptions$reltol, sqrt(.Machine$double.eps))
                            # steptol missing, 
                            # all parallel options missing
                  )
                 )
                 
                 
                 out$optim$phat <- out$optim$bestmem
                 out$optim$phattrans <- transFun(out$optim$bestmem, 
                                                 trans.L = returnL$transL$pretrans)
                 
                 
                 },
                 modMCMC = {
                           
                 sel <- which(parL$psel == 1)
                 
                 par.lo <- transFun(plo_bd, transL$ptrans)[sel]
                 par.up <- transFun(pup_bd, transL$ptrans)[sel]
                 p.ini <- par.lo + abs(par.lo-par.up)/2
                 
                 if(!"jump"%in%names(ALGoptions)) {ALGoptions$jump <- NULL }
                 if(!"prior"%in%names(ALGoptions)) {ALGoptions$prior <- NULL }
                 if(!"var0"%in%names(ALGoptions)) {ALGoptions$var0 <- NULL }
                 if(!"wvar0"%in%names(ALGoptions)) {ALGoptions$wvar0 <- NULL }
                 if(!"n0"%in%names(ALGoptions)) {ALGoptions$n0 <- NULL }
                 if(!"burninlength"%in%names(ALGoptions)) {ALGoptions$burninlength <- 0.3 }
                 if(!"covscale"%in%names(ALGoptions)) {ALGoptions$covscale <- 2.4^2/length(p.ini) }
                 if(!"ntrydr"%in%names(ALGoptions)) {ALGoptions$ntrydr <- 1 }
                 if(!"verbose"%in%names(ALGoptions)) {ALGoptions$verbose <- TRUE }
                 if(!"drscale"%in%names(ALGoptions)) {ALGoptions$drscale <- NULL }
                 
                 
                 outMCMC <- modMCMC(
                           f = resFun, 
                          # model specific
                          p = p.ini,
                          lower = par.lo,
                          upper = par.up,
                          shpmodel = shpmodel,
                          retdata = retdata, 
                          condata = condata, 
                          pretrans = returnL$transL$pretrans[sel], 
                          method = LikModel,
                          weight = weight, 
                          trim.query = trim.query,
                          ivap.query = ivap, 
                          parL = parL
                          # algorithm specific
                          , prior = ALGoptions$prior
                          , var0 = ALGoptions$var0 
                          , wvar0 = ALGoptions$wvar0
                          , n0 = ALGoptions$n0
                          , niter = ALGoptions$itermax
                          , outputlength = ALGoptions$itermax
                          , burninlength = ALGoptions$burninlength
                          , updatecov = ALGoptions$itermax/10
                          , covscale = ALGoptions$covscale
                          , ntrydr = 1
                          , drscale = ALGoptions$drscale
                          , verbose = TRUE
                 )
                           
                 out$optim$phat <- outMCMC$bestpar
                 out$optim$phattrans <- transFun(outMCMC$bestpar, 
                                                 returnL$transL$pretrans[sel] )
                 
                 },
                 stop("Enter a valid code for ALG")
          )
          
          # add the set up to
          
          settings <- vector("list")
          settings$data$retdata <- retdata
          settings$data$condata <- condata
          settings$transL <- transL
          settings$ALGoptions <- ALGoptions
          settings$weight <- weight
          settings$parL <- as.list(parL)
          settings$shpmodel
          settings$ivap <- ivap
          settings$hclip <- hclip
          settings$LikModel <- LikModel

          
          
          if(ALG == "DE"){
                    return(list("settings" = settings, "out" = out ))
          }
          
          if(ALG == "modMCMC"){
                    return(list("settings" = settings, "out" = outMCMC))
          }
          
}
