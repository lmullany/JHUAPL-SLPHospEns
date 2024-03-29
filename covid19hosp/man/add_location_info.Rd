% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hosp_ensemble_data_processing_functions.R
\name{add_location_info}
\alias{add_location_info}
\title{Pull in location information for the location codes in COVID hub processed data}
\usage{
add_location_info(df, location_variable = "location", keepvars = NULL)
}
\arguments{
\item{df}{data frame/data.table for which location information is required}

\item{location_variable}{character string that indicates the variable in df on which to match
the covid19hosp internal location dataset
a data table /data.frame that will be the source of the location information}

\item{keepvars}{an optional character vector of column names to retain after merge}
}
\value{
Returns a df merged with the location information
}
\description{
This function takes any data frame with a location code and a source for the lookup table. The
lookup can either be a path to a csv file, which the function then attempts to read in using
\code{data.table::fread}, or can be data.frame/data.table. Both must have a common variable on
which the merge will be based; the name of this common variable can indicated in "location_var"
which has default value "location"
}
\examples{
add_location_information(ensemble,"../covid19-forecast-hub/data-locations/locations.csv"), location_variable="location")
}
