\name{doy}
\alias{doy}
\title{
Day of year and days since a start month
}
\description{
Converts month names to the day of year of the first day of the month, and to the number of days from the first day of a start month.
Used for grazing and field application schedules.
}
\usage{
doy(month, begin = 'January')
}
\arguments{
  \item{month}{
    Character vector with English month names (e.g., \code{"March"}).
  }
  \item{begin}{
    Start month (English month name). Days are counted from the first day of this month, continuing into the next year if \code{month} comes before \code{begin}.
  }
}
\value{
  A list with elements \code{doy} (day of year of the first day of each month, non-leap year) and \code{day} (days from the first day of \code{begin}).
}
\examples{
doy(c('March', 'October'))
doy('February', begin = 'September')
}
