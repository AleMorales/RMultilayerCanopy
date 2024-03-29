% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/julia_interface.R
\name{interpolate_meteo}
\alias{interpolate_meteo}
\title{Calculate environmental conditions on a particular time of the day}
\usage{
interpolate_meteo(weather = DailyWeather(), time = weather$DL/2)
}
\arguments{
\item{weather}{Daily weather object as generated by \code{DailyWeather}}

\item{time}{Time of the day relative to sunrise (i.e., between 0 and weather$DL)}
}
\value{
A Julia object collecting the environmental variables at the specified time of the day:
\describe{
\item{Ib0}{Direct PAR that reaches the top of the canopy (mol/m2/s)}
\item{Id0}{Diffuse PAR that reaches the top of the canopy (mol/m2/s)}
\item{beta}{Solar elevation angle (rad)}
\item{Omega}{Solar azimuth angle (rad)}
\item{Tleaf}{Leaf temperature assumed to be equal to air (K)}
\item{VPD}{Vapor pressure deficit (Pa)}
\item{Ca}{CO2 molar fraction in the air (mol/mol)}
}
}
\description{
Calculate environmental conditions on a particular time of the day
}
\examples{
weather = DailyWeather()
env = interpolate_meteo(weather, weather$DL/2)
}
