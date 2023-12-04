\name{temperatures}
\alias{temperatures}
\alias{weather_dat_DK}
\alias{outside_slurry_temp_NIR}
\alias{outside_slurry_temp_dig_NIR}
\alias{outside_slurry_temp_vechi}
\alias{outside_slurry_temp_dig_vechi}

\docType{data}

\title{
Outside slurry temperatures and weather
}

\description{
Slurry temperatures according to the Danish National Inventory Report and from Vechi et al. 2023, and some related weather data for Denmark.
}

\usage{
data(weather_dat_DK)
data(outside_slurry_temp_NIR)
data(outside_slurry_temp_dig_NIR)
data(outside_slurry_temp_vechi)
data(outside_slurry_temp_dig_vechi)
}

\format{

These are temperature datasets for outside slurry storage in Denmark. 
The original data are from ...
\code{weather_dat_DK}: Danish weather data from 2000-2020 averaged by month.
\code{outside_slurry_temp_NIR}: Slurry temperatures by month in outside storages according to the Danish National Inventory Report 2022. The slurry temperature is estimated as 0.511 * air_temp + 5.1886.
\code{outside_slurry_temp_dig_NIT}: Slurry temperatures by month in outside storages receiving digestate according to the Danish National Inventory Report 2022. The slurry temperature is estimated as 0.75 * air_temp + 6.23.
\code{oustide_slurry_temp_vechi}: Measured slurry temperatures in pig slurry tanks in XXXXX and averaged by month. Data was averaged from Vechi et al. 2023. DOI: 10.1016/j.scitotenv.2023.166610.
\code{outside_slurry_temp_dig_vechi}: Measured slurry temperatures in slurry tanks recieving digestate in XXXX and averaged by month. Data was averaged from Vechi et al. 2023. DOI: 10.1016/j.scitotenv.2023.166610.
}

\source{

Nielsen, O.K., et al. 2023. Denmark's National Inventory Report 2023. DCE, Aarhus University. URL: https://dce2.au.dk/pub/SR541.pdf.

Vechi, N.T., Falk, J.M., Fredenslund, A.M., Edjabou, M.E., Scheutz, C., 2023. Methane emission rates averaged over a year from ten farm-scale manure storage tanks. Sci Total Environ. \doi{10.1016/j.scitotenv.2023.166610

}


}

\examples{
}

\keyword{datasets}
