% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intermittentATA-package.R, R/intermittent_ata.r
\docType{package}
\name{intermittentATA-package}
\alias{intermittentATA}
\alias{intermittentATA-package}
\title{intermittentATA: Intermittent demand time series analysis using the Ata Method based on the 'ATAforecasting package with the 'fable' Framework.}
\description{
Intermittent Ata Method is based on Croston's (1972) \doi{10.2307/3007885} method for intermittent demand forecasting, also described in Shenstone and Hyndman (2005) \doi{10.1002/for.963}. In Croston's approach, the non-zero parts of the time series are smoothed using simple exponential smoothing (SES), and the times between the non-zero items are smoothed separately. Allows ATA (Automatic Time series analysis using the Ata method) models from the 'ATAforecasting' package to be used in a tidy workflow with the modeling interface of 'fabletools'. This extends 'ATAforecasting' to provide enhanced model specification and management, performance evaluation methods, and model combination tools. The Ata method (Yapar et al. (2019) \doi{10.15672/hujms.461032}), an alternative to exponential smoothing (described in Yapar (2016) \doi{10.15672/HJMS.201614320580}, Yapar et al. (2017) \doi{10.15672/HJMS.2017.493}), is a new univariate time series forecasting method which provides innovative solutions to issues faced during the initialization and optimization stages of existing forecasting methods. Forecasting performance of the Ata method is superior to existing methods both in terms of easy implementation and accurate forecasting. It can be applied to non-seasonal or seasonal time series which can be decomposed into four components (remainder, level, trend and seasonal).

Intermittent Ata Method is based on Croston's (1972) \doi{10.2307/3007885} method for intermittent demand forecasting, also described in Shenstone and Hyndman (2005) \doi{10.1002/for.963}. In Croston's approach, the non-zero parts of the time series are smoothed using simple exponential smoothing (SES), and the times between the non-zero items are smoothed separately. Allows ATA (Automatic Time series analysis using the Ata method) models from the 'ATAforecasting' package to be used in a tidy workflow with the modeling interface of 'fabletools'. This extends 'ATAforecasting' to provide enhanced model specification and management, performance evaluation methods, and model combination tools. The Ata method (Yapar et al. (2019) \doi{10.15672/hujms.461032}), an alternative to exponential smoothing (described in Yapar (2016) \doi{10.15672/HJMS.201614320580}, Yapar et al. (2017) \doi{10.15672/HJMS.2017.493}), is a new univariate time series forecasting method which provides innovative solutions to issues faced during the initialization and optimization stages of existing forecasting methods. Forecasting performance of the Ata method is superior to existing methods both in terms of easy implementation and accurate forecasting. It can be applied to non-seasonal or seasonal time series which can be decomposed into four components (remainder, level, trend and seasonal).
}
\seealso{
Useful links:
\itemize{
  \item \url{https://alsabtay.github.io/intermittentATA/}
  \item Report bugs at \url{https://github.com/alsabtay/intermittentATA/issues}
}


Useful links:
\itemize{
  \item \url{https://alsabtay.github.io/intermittentATA/}
  \item Report bugs at \url{https://github.com/alsabtay/intermittentATA/issues}
}

}
\author{
\strong{Maintainer}: Ali Sabri Taylan \email{alisabritaylan@gmail.com} (\href{https://orcid.org/0000-0001-9514-934X}{ORCID}) [copyright holder]

Authors:
\itemize{
  \item Hanife Taylan Selamlar (\href{https://orcid.org/0000-0002-4091-884X}{ORCID}) [copyright holder]
  \item Guckan Yapar (\href{https://orcid.org/0000-0002-0971-6676}{ORCID}) [thesis advisor, copyright holder]
}

}
\keyword{internal}
\keyword{package}
