% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intermittent_ata.r
\name{glance.intermittentATA}
\alias{glance.intermittentATA}
\title{Glance an intermittentATA model}
\usage{
\method{glance}{intermittentATA}(x, ...)
}
\arguments{
\item{x}{An estimated model.}

\item{...}{Unused.}
}
\value{
A one row tibble summarising the model's fit.
}
\description{
Glance an intermittentATA model
}
\examples{
library(intermittentATA)
as_tsibble(fmcgData) \%>\%
   model(crostonata = intermittentATA(value ~ d_trend(type = "M", parQ = 1) + i_trend("A") + intermittent("croston"))) \%>\% glance()

}
