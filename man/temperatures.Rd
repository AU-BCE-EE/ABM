\name{temperatures}
\alias{temperatures}
\docType{data}
\title{
Outside slurry temperatures and weather
}

\description{
Slurry temperatures according to the Danish National Inventory Report and from Vechi et al. 2023.
}

\usage{
data(weather_dat_DK), Danish weather data from 2000-2020 averaged by month
data(outside_slurry_temp_NIR), slurry temperatures by month in outside storages according to the Danish National Inventory Report 2022. The slurry temperature is estimated as 0.511 * air_temp + 5.1886.
data(outside_slurry_temp_dig_NIR), slurry temperatures by month in outside storages receiving digestate according to the Danish National Inventory Report 2022. The slurry temperature is estimated as 0.75 * air_temp + 6.23.

data(outside_slurry_temp_vechi), Measured slurry temperatures in pig slurry tanks and averaged by month. Data was averaged from Vechi et al. 2023. DOI: 10.1016/j.scitotenv.2023.166610.
data(outside_slurry_temp_dig_vechi), Measured slurry temperatures in slurry tanks recieving digestate and averaged by month. Data was averaged from Vechi et al. 2023. DOI: 10.1016/j.scitotenv.2023.166610.

}

\format{
The format is: \rda

}
}
\source{

}

\examples{
plot(data(outside_slurry_temp_NIR))
viewSpec(outside_slurry_temp_vechi)
}

\keyword{datasets}