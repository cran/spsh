#' Calculates the Isothermal Water Vapour Conductivity
#' 
#' @description Calculates the isothermal vapour conductivity as a function of modelled volumetric air content. Different
#' models are implemented enabling the calculation of the relative gas diffusion coefficient (Ds/Do), based on 
#' different expressions for an effective tortuosity. 
#' 
#' @param p vector of soil hydraulic property model parameters, cf resp soil hydraulic property model for details.
#' @param por skalar value giving the fraction of a porous' media porosity [ - ]( value between [0, 1] ), defaults to the saturated water content. 
#' @param retFun soil hydraulic property function has to be specified if models \code{PMQ}, \code{TPM} or \code{TPEM} are used, necessary to calculate the air content at \code{h = 100} cm for the parameter \code{eps100}.
#' @param theta vector of numerical volumetric water contents [0,1] at which the air content is to be calculated.
#' @param model Implemented models (specify as character):
#' \tabular{lll}{
#'       \code{B}\tab{Buckingham (1904) }\cr
#'       \code{P}\tab{Penman (1940) }\cr
#'       \code{MQ60}\tab{ Millington and Quirck (1960) }\cr
#'       \code{MQ61}\tab{ Millington and Quirck (1961)}\cr
#'       \code{GS}\tab{Grable and Siemer (1968) }\cr
#'       \code{L}\tab{Lai et al. (1976) }\cr
#'       \code{PMQ}\tab{Moldrup et al. (1997) }\cr
#'       \code{TPM}\tab{Moldrup et al. (2004) }\cr
#'       \code{TPEM}\tab{Moldrup et al. (2005) }\cr
#'       
#' }
#' @param Temp Soil tempereature [ deg C ], defaults to 20.
#' @param m \code{PMQ} model parameter, default \code{m = 3}.
#' @param pF monotonically increasing pF values, defined as log10(| pressure head [ cm ]).
#' @param output Defaults to \code{log10} indicates the isothermal vapour conductivity is returned as \ifelse{html}{\out{log<sub>10</sub>}}{\eqn{log_10}}(conductivity), if \code{ouput != log10},
#' the output will be in non-transformed values.
#' @param ...  more arguments to be passed to \code{retFun}
#' 
#' @details More reading on the models reference is made to suggested in \insertCite{Weber.2019}{spsh}
#' 
#' @references 
#' \insertRef{Weber.2019}{spsh}\cr\cr
#' {{\bold{Buckingham, E.} (1904). Contributions to Our Knowledge of the Aeration Status of Soils, Bulletin 25, USDA Bureau of Soils, Washington, DC.}}
#' 
#' {{\bold{Grable, A.R.; Siemer, E.G.} (1968).Effects of Bulk Density, Aggregate Size, and Soil Water Suction on Oxygen Diffusion, Redox Potentials, and Elongation of Corn Roots. Soil Sci. Soc. Am. Proc., 32, pp. 180-186. <doi:10.2136/sssaj1968.03615995003200020011x>. }}
#' 
#' {{\bold{Lai, S.H.; Tiedje J.M.; Erickson, E.} (1976). In situ Measurement of Gas Diffusion Coefficient in Soils. Soil Sci. Soc. Am. J., 40, pp. 3-6. <doi:10.2136/sssaj1976.03615995004000010006x>. }}
#' 
#' {{\bold{Moldrup, P.; Olesen, T.; Rolston, D.E.; and Yamaguchi, T.} (1997). Modeling Diffusion and Reaction in Soils: Vii. Predicting Gas and Ion Diffusivity in Undisturbed and Sieved Soils. Soil Science. 162 (9): pp. 632-640. }}
#' 
#' {{\bold{Moldrup, P.; Olesen, T.; Yoshikawa, S.; Komatsu, T.; and Rolston, D.E.} (2004). Three-Porosity Model for Predicting the Gas Diffusion Coefficient in Undisturbed Soil. Soil Sci. Soc. Am. J. 68 (3).pp. 750-759. <doi:10.2136/sssaj2004.7500>.}}
#' 
#' {\bold{Moldrup, P.; Olesen, T.; Yoshikawa, S.; Komatsu, T.; and Rolston, D.E.} (2005). Predictive-Descriptive Models for Gas and Solute Diffusion Coefficients in Variably Saturated Porous Media Coupled to Pore-Size Distribution: II. Gas Diffusivity in Undisturbed Soil. Soil Sci., 170, pp. 854-866. <doi:10.1097/01.ss.0000196768.44165.1f>.}
#' 
#' {{\bold{Millington, R.J.; Quirk, J.P.} (1960). Millington, R. J., and Quirk. J.M.
#'       Transport in porous media. pp. 97-106. In: F.A. Van Beren, et al. (ed.) Trans. Int. Congr. Soil Sci., 7 th,
#'       Vol. 1, Madison, Wl. 14-24 Aug. 1960. Elsevier, Amsterdam.}}
#' 
#' {{\bold{Millington, R.J.; Quirk, J.P.} (1961). Permeability of Porous Solids. Trans. Faraday Soc., 1961, 57, pp. 1200-1207. <doi:10.1039/TF9615701200>.}}
#' 
#' {{\bold{Penman, H.L.} (1940). Gas and vapour movements in the soil: I. The diffusion of vapours through porous solids. J. Agric. Sci., 30: pp. 437-462. <doi:10.1017/S0021859600048164>. }}
#' 
#' 
#' {{\bold{Xu, X; Nieber, J.L. Gupta, S.C.} (1992). Compaction Effect on the Gas Diffusion Coefficient in Soils. Soil Sci. Soc. Am. J.,56, pp. 1743-1750. <doi:10.2136/sssaj1992.03615995005600060014x>.}}
#' 
#' 
#' @author Tobias KD Weber , \email{tobias.weber@uni-hohenheim.de}
#'
#' @examples
#' # | pressure head |
#' pF <- seq(-3, 7, length = 201)
#' h <- 10^pF
#' # van Genuchten-Mualem model parameters
#' p <- c(0.08, .42, .05, 1.5, 100, .5)
#' # calculate soil hydraulic property values
#' shypL <- shypFun.01110(p, h)
#' # clculate the isothermal vapour conductivity
#' kvap <- KvapFun(p, por = p[2], retFun = NA, theta = shypL$theta, model = "MQ61", 
#'                 Temp = 20, m = 3, pF, output = "log10")
#' 
#' @importFrom Rdpack reprompt

#' @export

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
