

#' Specify parameters of the virtual crop
#'
#' @param Sco25 Rubisco specificity factor at 25 C
#' @param E_Sco Apparent activation energy of Sc/o (kJ/mol)
#' @param Kmc25 Michaelis-Menten constant of carboxylation with respect to CO2 at 25 C (mol/mol)
#' @param E_Kmc Activation energy of Kmc (kJ/mol)
#' @param Kmo25 Michaelis-Menten constant of oxygenation with respect to O2 at 25 C (mol/mol)
#' @param E_Kmo Activation energy of Kmo (kJ/mol)
#' @param ar Ratio between Vcmax25 and Nr (1/s)
#' @param E_Vcmax Activation energy of Vcmax (kJ/mol)
#' @param theta Curvature parameter in the light response of the electron transport
#' @param Knc Nc at which the leaf absorbs 50\% of incoming PAR (mol/m2)
#' @param Phi2LL PSII quantum yield at low light
#' @param Phi1LL PSI quantum yield at low light
#' @param fcyc Fraction of cyclic electron transport
#' @param Topt_Phi2 Optimal temperature for PSII quantum yield
#' @param Omega sigma/√2 in the Gaussian temperature function for PSII quantum yield
#' @param aj Ratio between Jmax25 and Nt (1/s)
#' @param ks Ratio between Ns and Jmax25 (s)
#' @param E_Jmax Activation energy Jmax (kJ/mol)
#' @param D_Jmax Deactivation energy of Jmax (kJ/mol)
#' @param S_Jmax Entropy coefficient of Jmax (kJ/mol/K)
#' @param E_Rd Activation energy of Rd (kJ/mol)
#' @param f_Rd Ratio between Rd25 and Vcmax25
#' @param gm25 Mesophyll conductance (mol/m2/s)
#' @param E_gm Activation energy gm (kJ/mol)
#' @param S_gm Entropy coefficient of gm (kJ/mol/K)
#' @param D_gm Deactivation energy of gm (kJ/mol)
#' @param gs0 Minimum stomatal conductance (mol/m2/s)
#' @param sgs Scaling factor between gross assimilation and stomatal conductance
#' @param D0  Sensitivity of gs to VPD (kPa)
#' @param gb Boundary layer conductance (mol/m2/s)
#' @param angles Named list with parameters for the leaf angle distribution (see documentation)
#' @param Ncmin Minimum leaf nitrogen content (g N/g DW)
#' @param sigma_soil Albedo of the soil surface
#' @param sky List with parameters for the sky model (see documentation)
#' @param lat Latitude of the location in degrees (positive values for Northern hemisphere)
#' @param Ca Air CO2 molar ratio (umol/mol)
#'
#' @details
#' The argument \code{angles} should be a named list with components \code{X} and\code{rtol}, the
#' first one being the parameter of the ellipsoidal leaf angle distribution model (a positive value)
#' and the second being the relative tolerance used by the adaptive cubature algorithm used to
#' calculated CO2 assimilation in the sunlit leaves.
#'
#' The argument \code{sky} should be a named list with components \code{mode} and \code{nsky},
#' the first one being the string \code{"standard"} or \code{"uniform"} depending on the type of
#' sky model to be used and \code{nsky} being the number of hemispherical rings into which the
#' sky should be divided.
#'
#' @return
#' A Julia object with all the parameter values (access via $ notation)
#'
#' @export
#'
#' @examples
#' pars = parameters()
#' pars$lat
parameters = function(
  Sco25 = 2800.0,
  E_Sco = -24.46,
  Kmc25 = 270.0,
  E_Kmc = 80.99,
  Kmo25 = 165.0,
  E_Kmo = 23.72,
  ar =3.5/786,
  E_Vcmax = 65.33,
  theta = 0.7,
  Knc = 0.076/25.0,
  Phi2LL = 0.85,
  Phi1LL = 1.0,
  fcyc = 0.1,
  Topt_Phi2 = 22.5 + 273.15,
  Omega = 36.5,
  aj = 15.9e-3,
  ks = 125.0,
  E_Jmax = 30.0,
  D_Jmax = 200.0,
  S_Jmax = 0.65,
  E_Rd = 46.39,
  f_Rd  = 0.01,
  gm25 = 0.4,
  E_gm = 70.2,
  S_gm = 0.32,
  D_gm = 94.0,
  gs0 = 0.05/1.56,
  sgs = 3.0,
  D0 = 1.5,
  gb = 0.5,
  angles = list(X = 1.0, rtol = 1e-3),
  Ncmin   = 8.4e-3,
  sigma_soil = 0.21,
  sky = list(mode = "standard", nsky = 10),
  lat = 52.0,
  Ca = 4e2) {

  JuliaCall::julia_library("MultilayerCanopy")
  JuliaCall::julia_eval("const MC = MultilayerCanopy")

  # Create LeafAngleModel
  lambda = JuliaCall::julia_call("MC.Ellipsoidal", angles$X, need_return = "Julia")
  leafangle = JuliaCall::julia_call("MC.LeafAngleModel", lambda, angles$rtol)

  # Create modified sky model
  sky_model = ifelse(sky$mode == "standard", "MC.standard_sky", "MC.uniform_sky")
  sky_n = JuliaCall::julia_call("Val", as.integer(sky$nsky), need_return = "Julia")
  rawsky = JuliaCall::julia_call(sky_model, sky_n, need_return = "Julia")
  upsky = JuliaCall::julia_call("MC.ksky", leafangle, rawsky, need_return = "Julia")

  JuliaCall::julia_call("MC.parameters",
                        Sco25 = Sco25,
                        E_Sco = E_Sco,
                        Kmc25 = Kmc25,
                        E_Kmc = E_Kmc,
                        Kmo25 = Kmo25,
                        E_Kmo = E_Kmo,
                        ar = ar,
                        E_Vcmax = E_Vcmax,
                        theta = theta,
                        Knc = Knc,
                        Phi2LL = Phi2LL,
                        Phi1LL = Phi1LL,
                        fcyc = fcyc,
                        Topt_Phi2 = Topt_Phi2,
                        Omega = Omega,
                        aj = aj,
                        ks = ks,
                        E_Jmax = E_Jmax,
                        D_Jmax = D_Jmax,
                        S_Jmax = S_Jmax,
                        E_Rd = E_Rd,
                        f_Rd  = f_Rd,
                        gm25 = gm25,
                        E_gm = E_gm,
                        S_gm = S_gm,
                        D_gm = D_gm,
                        gs0 = gs0,
                        sgs = sgs,
                        D0 = D0,
                        gb = gb,
                        angles = leafangle,
                        Ncmin   = Ncmin,
                        sigma_soil = sigma_soil,
                        sky = upsky,
                        lat = lat,
                        Ca = Ca)
}



#' Specify the variables of the virtual crop
#'
#' @param pars Parameters as generated by the \code{parameters}
#' @param nleaf Canopy nitrogen content (g N/m2)
#' @param wleaf Canopy (leaf) biomass as dry weight (g/m2)
#' @param SLA Specific leaf area (cm2/g)
#' @param kN Extinction coefficient of leaf nitrogen within the canopy
#' @param f_Ncmn Minimum fraction of N allocated to chlorophyll
#' @param f_Ncmx Maximum fraction of N allocated to chlorophyll
#' @param pf_Nc Scaling constant for fraction of N allocated to chlorophyll
#' @param f_Nrmn Minimum fraction of N allocated to Rubisco
#' @param f_Nrmx Maximum fraction of N allocated to Rubisco
#' @param pf_Nr Scaling constant for fraction of N allocated to Rubisco
#'
#' @return
#' A Julia object with all the variables required to perform the simulations:
#' \describe{
#' \item{Nmin}{Minimum leaf nitrogen content (mol/m2)}
#' \item{GLAI}{Green leaf area index (limited by \code{nleaf} or \code{wleaf})}
#' \item{Nlmax}{Leaf nitrogen content on top of the canopy (mol/m2)}
#' \item{kN}{Extinction coefficient of leaf nitrogen within the canopy}
#' \item{f_Ncmn}{Minimum fraction of N allocated to chlorophyll}
#' \item{f_Ncmx}{Maximum fraction of N allocated to chlorophyll}
#' \item{pf_Nc}{Scaling constant for fraction of N allocated to chlorophyll}
#' \item{f_Nrmn}{Minimum fraction of N allocated to Rubisco}
#' \item{f_Nrmx}{Maximum fraction of N allocated to Rubisco}
#' \item{pf_Nr}{Scaling constant for fraction of N allocated to Rubisco}
#' }
#'
#' @export
#'
#' @examples
#' pars = parameters()
#' vars = variables(pars)
#' vars$GLAI
variables = function(
                      pars,
                      nleaf  = 4.0,
                      wleaf = 100.0,
                      SLA = 300.0,
                      kN   = 0.4,
                      f_Ncmn = 0.1,
                      f_Ncmx = 0.3,
                      pf_Nc  = 4.0,
                      f_Nrmn = 0.2,
                      f_Nrmx = 0.6,
                      pf_Nr  = 2.0) {

  JuliaCall::julia_library("MultilayerCanopy")
  JuliaCall::julia_eval("const MC = MultilayerCanopy")
  JuliaCall::julia_call("MC.variables",
                        pars,
                        nleaf  = nleaf,
                        wleaf = wleaf,
                        SLA = SLA,
                        kN = kN,
                        f_Ncmn = f_Ncmn,
                        f_Ncmx = f_Ncmx,
                        pf_Nc  = pf_Nc,
                        f_Nrmn = f_Nrmn,
                        f_Nrmx = f_Nrmx,
                        pf_Nr  = pf_Nr,
                        need_return = "Julia")

}


#' Convert hexadecimal degrees to radians
#'
#' @param x Angle in hexadecimal degrees
#'
#' @return
#' Angle in radians
#' @export
#'
#' @examples
#' toRad(52)
toRad = function(x) {
  x*pi/180
}

#' Combine weather information at the daily level
#'
#' @param lat Latitude of the location in hexadecimal degrees
#' @param doy Day of the year being simulated
#' @param tau Atmospheric transmissivity
#' @param Tmin Minimum temperature of the day (Celsius)
#' @param Tmax Maximum temperature of the day (Celsius)
#' @param RH Relative humidity at sunrise
#' @param Ca CO2 molar fraction in the air (umol/mol)
#'
#' @return
#' A Julia object collecting the weather variables after some transformations:
#' \describe{
#' \item{lat}{Latitude of the location (rad)}
#' \item{doy}{Day of the year being simulated}
#' \item{tsr}{Time of the day at which sunrise occurs (h)}
#' \item{DL}{Diurnal length, that is, difference between sunrise and sunset (h)}
#' \item{tau}{Atmospheric transmissivity}
#' \item{Temp}{Function to interpolate temperature during the day (does not work in R)}
#' \item{ea}{Air vapour pressure at sunrise (kPa)}
#' \item{Ca}{CO2 molar fraction in the air (umol/mol)}
#' }
#' @export
#'
#' @examples
#' weather = DailyWeather()
#' weather$DL
DailyWeather = function(lat = 52, doy = 182L, tau = 0.75, Tmin = 20,
                        Tmax = 30, RH = 0.4, Ca = 400) {
  JuliaCall::julia_library("MultilayerCanopy")
  JuliaCall::julia_eval("const MC = MultilayerCanopy")
  JuliaCall::julia_call("MC.DailyWeather", toRad(lat), doy, tau, Tmin, Tmax, RH, Ca,
                        need_return = "Julia")
}


#' Calculate environmental conditions on a particular time of the day
#'
#' @param weather Daily weather object as generated by \code{DailyWeather}
#' @param time Time of the day relative to sunrise (i.e., between 0 and weather$DL)
#'
#' @return
#' A Julia object collecting the environmental variables at the specified time of the day:
#' \describe{
#' \item{Ib0}{Direct PAR that reaches the top of the canopy (mol/m2/s)}
#' \item{Id0}{Diffuse PAR that reaches the top of the canopy (mol/m2/s)}
#' \item{beta}{Solar elevation angle (rad)}
#' \item{Omega}{Solar azimuth angle (rad)}
#' \item{Tleaf}{Leaf temperature assumed to be equal to air (K)}
#' \item{VPD}{Vapor pressure deficit (Pa)}
#' \item{Ca}{CO2 molar fraction in the air (mol/mol)}
#' }
#' @export
#'
#' @examples
#' weather = DailyWeather()
#' env = interpolate_meteo(weather, weather$DL/2)
interpolate_meteo = function(weather = DailyWeather(), time = weather$DL/2) {
  JuliaCall::julia_library("MultilayerCanopy")
  JuliaCall::julia_eval("const MC = MultilayerCanopy")
  JuliaCall::julia_call("MC.interpolate_meteo", weather, time,
                        need_return = "Julia")
}


#' Compute CO2 assimilation for a leaf given parameters and environment
#'
#' @param pars Parameters as generated by \code{\link{parameters}}
#' @param Np Photosynthetic nitrogen content (mol/m2)
#' @param f_Nc Fraction of photosynthetic nitrogen allocated to chlorophyll
#' @param f_Nr Fraction of photosynthetic nitrogen allocated to Rubisco
#' @param VPD Vapor pressure deficit (kPa)
#' @param PAR Photosynthetically active radiation (umol/m2/s)
#' @param Tleaf Leaf temperature (Celsius)
#' @param Ca CO2 molar fraction in the air (umol/mol)
#'
#' @return
#' A named vector with gross CO2 assimilation (\code{"Ag"}), net CO2 assimilation (\code{"A"}),
#' net CO2 assimilation limited by carboxylation capacity (\code{"Ac"}) and
#' net CO2 assimilation limited by electron transport capacity (\code{"Aj"})
#' @export
#'
#' @examples
#' pars = parameters()
#' As = A(pars)
#' As["A"] == min(As["Ac"], As["Aj"])
A = function(pars = parameters(), Np = 0.1, f_Nc = 0.15, f_Nr = 0.30,
              VPD = 1.0, PAR = 1e3, Tleaf = 25, Ca = 4e2) {
  JuliaCall::julia_library("MultilayerCanopy")
  JuliaCall::julia_eval("const MC = MultilayerCanopy")
  out = JuliaCall::julia_call("MC.Ag", pars, Np, f_Nc, f_Nr, VPD*1e3,
                              PAR/1e6, Tleaf + 273.15, Ca/1e6,
                              need_return = "R")
  out = as.numeric(out)*1e6
  names(out) = c("Ag", "A", "Ac", "Aj")
  out
}


#' Ellipsoidal leaf angle distribution
#'
#' @param X Ellipsoidal parameters (must have positive value)
#' @param angles Leaf inclination angles (radians)
#'
#' @return
#' Probability density of the leaf inclination angle(s)
#' @export
#'
#' @examples
#' curve(ellipsoidal(2.5, x), 0, pi/2)
ellipsoidal = function(X, angles) {
  if(X < 1) {
    eps = sqrt(1 - X^2)
    Lambda = X + asin(eps)/eps
  } else if(X == 1) {
    Lambda = 2
  } else {
    eps = sqrt(1 - 1/X^2)
    Lambda = X + log((1 + eps)/(1 - eps))/(2*eps*X)
  }
  den = (cos(angle)^2 + X^2*sin(angle)^2)
  2*X^3*sin(angle)/(Lambda*den^2)
}


#' Canopy structure in layers
#'
#' @param pars Parameters as generated by \code{\link{parameters}}
#' @param vars Variables as generated by \code{\link{variables}}
#' @param nlayers Number of layers into which the canopy should be divided
#'
#' @return
#' A Julia object with photosyntetic traits per layer as well as pre-calculated
#' values required for computing PAR distribution
#' \describe{
#' \item{DeltaL}{Increment of leaf area index of each layer}
#' \item{Lm}{Cumulative LAI from the top of the canopy to the middle of each layer}
#' \item{Lu}{Cumulative LAI from the top of the canopy to the upper boundary of each layer}
#' \item{Fu}{Matrix to compute fraction of upward scattered PAR that is intercepted}
#' \item{Fd}{Matrix to compute fraction of downward scattered PAR that is intercepted}
#' \item{fdif}{Fraction of diffuse PAR intercepted by each layer}
#' \item{Vcmax25}{Temperature-normalized maximum rate of carboxylation in each layer}
#' \item{Rd25}{Temperature-normalized of mitochondrial respiration in each layer}
#' \item{Jmax25}{Temperature-normalized maximum rate of electron transport in each layer}
#' \item{alpha}{Leaf absorptance in each layer}
#' \item{sigma}{Leaf scattering coefficient in each layer}
#' }
#' @export
#'
#' @examples
#' pars = parameters()
#' vars = variables(pars)
#' can = canopy(pars, vars, 15)
canopy = function(pars, vars, nlayers) {
  JuliaCall::julia_library("MultilayerCanopy")
  JuliaCall::julia_eval("const MC = MultilayerCanopy")
  nlay = JuliaCall::julia_call("Val", as.integer(nlayers),need_return = "Julia")
  can = JuliaCall::julia_call("MC.Canopy", pars, vars, nlay,
                              need_return = "Julia")
}

#' Calculate canopy photosynthesis and intermediate variables per layer
#'
#' @param can Canopy structure object as generated by \code{\link{canopy}}
#' @param env Object containing environmental conditions at a particular time as
#' generated by \code{\link{interpolate_meteo}}
#' @param pars Parameters as generated by \code{\link{parameters}}
#' @param vars Variables as generated by \code{\link{variables}}
#'
#' @return
#' List with the following information:
#' \describe{
#' \item{Acan}{Total canopy gross CO2 assimilation per unit of ground (mumol/m2/s)}
#' \item{Alayer}{Canopy gross CO2 assimilation of each layer per unit of ground (mumol/m2/s)}
#' \item{Ilayer}{Average PAR intercepted by each layer per unit of ground (mumol/m2/s)}
#' \item{fsun}{Fraction of leaves that are sunlit in each layer}
#' \item{Ashade}{Average net CO2 assimilation per unit of shaded leaf area in each layer (mumol/m2/s)}
#' \item{Asun}{Average net CO2 assimilation per unit of sunlit leaf area in each layer (mumol/m2/s)}
#' }
#' @export
#'
#' @examples
#' pars = parameters()
#' vars = variables(pars)
#' can = canopy(pars, vars, 15)
#' weather = DailyWeather()
#' env = interpolate_meteo(weather, weather.DL/2)
#' can_fluxes = Acan(can, env, pars, vars)
#' plot(can$Lm, can_fluxes$Alayer)
Acan = function(can, env, pars, vars) {
  JuliaCall::julia_library("MultilayerCanopy")
  JuliaCall::julia_eval("const MC = MultilayerCanopy")
  out = JuliaCall::julia_call("MC.Acan", can, env, pars, vars, need_return = "R")
  names(out)= c("Acan", "Alayer", "Ilayer", "fsun", "Ashade", "Asun")
  out$Acan = out$Acan*1e6
  out$Alayer = out$Alayer*1e6
  out$Ilayer = out$Ilayer*1e6
  out$Ashade = out$Ashade*1e6
  out$Asun = out$Asun*1e6
  out
}


#' Calculate diurnal net CO2 assimilation
#'
#' @param weather Object with daily weather variables as generated by \code{\link{DailyWeather}}
#' @param can Canopy structure object as generated by \code{\link{canopy}}
#' @param pars Parameters as generated by \code{\link{parameters}}
#' @param vars Variables as generated by \code{\link{variables}}
#' @param nsteps Numbers of segments in which the day should be subdivided for numerical integration
#' @param ngrid Order of the Gaussian-Legendre integration to be performed within each segment
#'
#' @return
#' @export
#'
#' @examples
#' weather = DailyWeather()
#' pars = parameters()
#' vars = variables(pars)
#' can = canopy(pars, vars, 15)
#' Aday = calc_Adiurnal(weather, can, pars, vars)
calc_Adiurnal = function(weather, can, pars, vars, nsteps = 3, ngrid = 5) {
  JuliaCall::julia_library("MultilayerCanopy")
  JuliaCall::julia_eval("const MC = MultilayerCanopy")
  ngrid = JuliaCall::julia_call("Val", as.integer(ngrid))
  nsteps = JuliaCall::julia_call("Val", as.integer(nsteps))
  method = JuliaCall::julia_call("MC.PGLIntegrator", ngrid, nsteps)
  JuliaCall::julia_call("MC.calc_Adiurnal", weather, can, pars, vars, method,
                        need_return = "R")
}



