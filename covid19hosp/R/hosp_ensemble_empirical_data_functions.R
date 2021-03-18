#' Pull empirical hospitalization data
#'
#' As of 2020/11/16, healthdata.gov began providing publicly available
#' csv formatted data of hospitalizations associated with covid. This
#' function takes a url, reads in the timeseries data, and converts it
#' to long format,combining pediatric and adult, and creating a US version
#' @param source string either "state" or "facility"; used internally by the function to pull the fixed
#' identifier for the dataset.
#'  - If "facility": not yet implemented
#'  - If "state" (default), the dataset at
#'  "https://healthdata.gov/api/views/g62h-syeh/rows.csv?accessType=DOWNLOAD&api_foundry=true"
#'  will be returned
#' @return data.table of hospitalization data
#' @export
#' @examples
#' pull_empirical_hospitalization_data()
#' pull_empirical_hospitalization_data(url)
#' pull_empirical_hospitalization_data(url=hosp_url)

pull_empirical_hospitalization_data <- function(source=c("state")) {

  source = match.arg(source)
  state_url = "https://healthdata.gov/api/views/g62h-syeh/rows.csv?accessType=DOWNLOAD&api_foundry=true"


  if(source=="state") raw <- pull_and_process_state(state_url)
  if(source=="facility") raw <- pull_and_process_facility(url)

  #ensure that date is data variable not character)
  raw[,date:=as.IDate(date)]

  #keep columns of interest, and set names
  empir <- raw[,.(state,date,confirmed_covid_24h,suspected_covid_24h,total_covid_24h)]
  beds <- raw[,.(state,date,inpatient_beds, inpatient_beds_used)]

  #melt this long
  empir <- melt(empir,id.vars = c("state","date"), variable.name="indicator",value.name="hospitalizations")

  #add the US to the empir
  empir <- rbind(empir, empir[,.("hospitalizations"=sum(hospitalizations,na.rm=T)),by=.(date,indicator)][,state:="US"])

  #add the US to the beds
  beds <- rbind(beds, beds[,.("inpatient_beds"=sum(inpatient_beds, na.rm=T),
                                "inpatient_beds_used" = sum(inpatient_beds_used, na.rm=T)),
                             by=.(date)][,state:="US"])
  # Change the covid_count names
  empirical_outcome_names = list("Confirmed" = "confirmed_covid_24h",
                                 "Suspected" = "suspected_covid_24h",
                                 "Confirmed + Suspected" = "total_covid_24h")

  empir[,indicator:=factor(indicator,
                          levels=empirical_outcome_names,
                          labels=names(empirical_outcome_names))]

  #Add the internal location dataset
  empir <- location_dataset[,.(state, state_name)][empir,on="state"]

  #subtract the date by one
  empir[, date:=date-1]

  return(merge(empir, beds, by=c("state","date"), all.x=T))


}


#' Construct hosp url
#'
#' Construct a time series url for pulling the latest time series data
#' @param source string either "state" or "facility"; used internally by the function to pull the fixed
#' identifier for the dataset. Only "state" is currently implemented.
#' @return string url to data; this can be fed to `readr::read_csv`, `datatable::fread()` or other similar function
#' @export
#' @examples
#' construct_hosp_url("state")
construct_hosp_url <- function(source = c("state","facility")) {

  source = match.arg(source)
  if(source=="facility") {
    stop("Not implemented yet for facility level data", call. = F)
  }

  id_list = list(
    #by facility dataset: https://healthdata.gov/dataset/covid-19-reported-patient-impact-and-hospital-capacity-facility
    "facility" = "d475cc4e-83cd-4c16-be57-9105f300e0bc",
    #by state dataset: https://healthdata.gov/dataset/covid-19-reported-patient-impact-and-hospital-capacity-state-timeseries
    "state" = "83b4a668-9321-4d8c-bc4f-2bef66c49050"
  )
  id = id_list[[source]]

  tsurl <- jsonlite::fromJSON(paste0("https://healthdata.gov/api/3/action/package_show?id=",id,"&page=0"))$result$resources[[1]]$url

  return(tsurl)
}


#' Function to pull the and pre-process the raw state level file
#'
#' Function will pull and process the raw state level file.
#' @param url The url for the data
#' @keywords internal
#'
pull_and_process_state <- function(url) {
  raw = try(data.table::fread(url),silent=T)
  if(class(raw) == "try-error") {
    cat("\ncurl download error / trying (slower) read.csv\n")
    raw = data.table::setDT(read.csv(url))
  }
  raw <- raw[,.SD,.SDcols = patterns("state|date|^previous.+(confirmed|suspected)$|inpatient_beds|inpatient_beds_used")]
  raw[,confirmed_covid_24h := rowSums(.SD,na.rm=T), .SDcols = c("previous_day_admission_adult_covid_confirmed","previous_day_admission_pediatric_covid_confirmed")]
  raw[,suspected_covid_24h := rowSums(.SD,na.rm=T), .SDcols = c("previous_day_admission_adult_covid_suspected","previous_day_admission_pediatric_covid_suspected")]
  raw[,total_covid_24h:=rowSums(.SD,na.rm=T),.SDcols = c("confirmed_covid_24h","suspected_covid_24h")]
  return(raw)
}

#' Function to pull the and pre-process the raw facility level weekly file
#'
#' Function will pull and process the raw facility level weekly file.
#' @param url The url for the data
#' @keywords internal
#'
pull_and_process_facility <- function(url) {
  #Not yet implemented
  return(NULL)
}

