# Ensemble Pre-processing function

#' Get raw data for a modeling group
#'
#' A function to get the last dataset, dated up to a given date, for a given modeling group. The function
#' will attempt to find the single matching csv file for this modeling group that is up-to, but
#' not exceeding the date provided
#' @param modelname character name of the model (must be subfolder in context)
#' @param date yyyy-mm-dd character indicating the date threshold for the model
#' @param context character path to the data-processed folder for covid19-forecast-hub
#' @return data table containing the data for this modeling group
#' @import data.table
#' @export
#' @examples
#' get_data_for_model(modelname="JHU-APL",date="2020-11-01")
#' get_data_for_model(modelname="JHU-APL",date="2020-11-01",context =mypath_to_covid_hub)
get_data_for_model <- function(modelname, date, context="../covid19-forecast-hub/data-processed/") {
  fnames <- dir(path = paste0(context,modelname),pattern = ".csv")

  if(length(fnames)<1) return(NULL)
  fdates <- as.Date(stringr::str_extract(fnames,"[0-9]{4}-[01][0-9]-[0123][0-9]"))
  if(length(which(fdates<=date))==0) return(NULL)
  maxdate <- fdates[max(which(fdates<=date))]
  lastfile <- fnames[stringr::str_detect(fnames,pattern=as.character(maxdate))]
  if(length(lastfile)!=1) {
    stop("More than one data file found?")
  } else {
    fullname <- paste0(context,modelname,"/",lastfile)
    DT = data.table::fread(fullname)[,src:=modelname]
    return(DT)
  }
}

#' Return a list of model directories
#'
#' Helper function to get the list of model groups in the covid19-hub processed data collection
#' @param context character path to the data-processed folder for covid19-forecast-hub
#' @return list of subfolders in the covid19 process data subfolder
#' @export
#' @examples
#' get_model_directories()
#' get_model_directories(context =mypath_to_covid_hub)
get_model_directories <- function(context="../covid19-forecast-hub/data-processed/") {
  model_dirs <- list.dirs(context, full.names=FALSE)[-1]
  # always exclude UMass ensembles and our own ensemble
  model_dirs <- model_dirs[!model_dirs %in% c(
    "JHUAPL-SLPHospEns",
    "COVIDhub-ensemble",
    "COVIDhub-baseline",
    "COVIDhub-trained_ensemble"
    )]
  return(model_dirs)
}

#' Get the hospitalization quantile functions by group
#'
#' This function is the workhorse of the pre-processing set of functions. Given a list of modeling
#' group subfolders, and a date threshold, it will pull the raw data (if available),
#' limit the data to hospitalization outcomes, and limit the dates to certain parameters
#' @param model_dir_list a list of character strings, one for each subfolder in the data-processed folder for covid19-forecast-hub
#' @param date_threshold a date threshold - only data prior to this date will be pulled; this data will also be used to restrict
#' the forecast date for any data sets returned. Only if the forecast date is greater than or equal to
#' six days prior to the forecast date, will the data be retained.
#' @param date_limits a two Date vector of limits on the data to be returned. The default is to return data from the threshold date
#' to 28 days in the future
#' @param context character path to the data-processed folder for covid19-forecast-hub
#' @param includepoint default FALSE (normally, we will drop the point estimate)
#' @return raw hospitalization quantile functions by modeling group
#' @export
#' @examples
#' get_hosp_raw_data(mydirs, "2020-11-02")
get_hosp_raw_data <- function(model_dir_list,
                              date_threshold,
                              date_limits=c(as.Date(date_threshold),as.Date(date_threshold)+28),
                              context="../covid19-forecast-hub/data-processed/",
                              includepoint = F) {

  # Pull ALL the data from all the modeling groups
  hosp <- data.table::rbindlist(lapply(model_dir_list,
                                       get_data_for_model,
                                       date=date_threshold,
                                       context=context), use.names = TRUE)

  #Make sure the location variable for states is two digits
  hosp[, location:=data.table::fifelse(stringr::str_length(location)<=2, stringr::str_pad(location,width=2, pad='0'),location)]

  #Update the target variable so that days ahead less than 10 get a leading zero
  hosp[stringr::str_detect(target,"^\\d "),target:=paste0("0",target)]

  # limit to hospital outcome
  hosp <- hosp[
    stringr::str_detect(target,"hosp") &
      forecast_date>=(as.Date(date_threshold)-6) &
      between(target_end_date, date_limits[1],date_limits[2])]

  # Normally, we will be excluding the point estimate, (i.e. restricting to "quantile")
  if(!includepoint) {
    hosp <- hosp[type=="quantile"]
  }

  data.table::setorder(hosp,src,location,target_end_date,quantile)

  return(hosp)

}

#' Check consistency of hospitalization data
#'
#' Function provides some checks on the hospitalization data
#' @param hosp Data table containing the modeling groups raw quantile functions by location and date
#' @param datethreshold string date "YYYY-MM-DD" indicating the forecast date
#' @param valid_models a list of valid hospitalization modeling groups. The default is to pull from
#' `valid_hospitalization_models()` which is a hard coded list available in the package
#' @param cu_select_only logical default = T include only the CU-select model if there are multiple Columbia Models
#' @param use_target_label_28 logical default = F, use the target label to isolate 28 days of forecast
#' @return TRUE if all tests are passed; otherwise error condition
#' @export
#' @examples
#' hosp_processing(hosp)
hosp_processing <- function(hosp, datethreshold, valid_models=valid_hospitalization_models(),cu_select_only=TRUE, use_target_label_28=FALSE) {

  if(cu_select_only) {
    hosp <- hosp[!src %chin% c("CU-nochange","CU-scenario_high", "CU-scenario_low", "CU-scenario_mid")]
  }

  # Check for any src values that are not in hosp_models
  if(length(setdiff(unique(hosp$src), valid_models))!=0) {
    print("There is model in the src column of hosp frame, that is not in the valid hospitalization list")
    stop(call. = F)
  }

  # Before checking for any varying dates, if use_target_label_28=T, then we will simply take all the 01
  # through 28 days ahead rows, and set the target_end date
  if(use_target_label_28) {
    mindate = hosp[,min(target_end_date)]
    hosp[,target_end_date:=mindate+as.integer(stringr::str_sub(target,1,2))-1]
    hosp <- hosp[between(target_end_date,as.Date(mindate), as.Date(mindate+27))]
  } else {
    # we will in this case, take the first 28 days available after the date threshold for each model
    hosp <- hosp[between(target_end_date,as.Date(datethreshold)+1, as.Date(datethreshold)+28)]
  }

  # Check for any varying dates
  if(!estimates_consistent(hosp)) {
    print("There are one or more locations that have varying model estimates by target_end_date")
    stop(call. = F)
  }

  # Check that quantile is numeric/double
  if(suppressWarnings(all(!is.na(as.numeric(hosp$quantile))))) {
    hosp[,quantile:=as.numeric(quantile)]

  } else {
      print("Quantile is not numeric, and cannot be converted")
    stop(call. = F)
  }

  # At this point, we can return the modified frame, since all tests passed
  return(hosp)
}


#' Get a current list of valid hospitalization models
#'
#' Function will return a list of current valid hospitalization models
#' @return list of valid hospitalization models
#' @export
#' @examples
#' valid_hospitalization_models()
valid_hospitalization_models <- function() {
  return(list("Columbia" = "CU-select",
              "Covid19Sim" = "Covid19Sim-Simulator",
              "ERDC" = "USACE-ERDC_SEIR",
              "GoogleHarvard" = "Google_Harvard-CPF",
              "GT-DeepCOVID" = "GT-DeepCOVID",
              "IHME" = "IHME-CurveFit",
              "JHUAPL" = "JHUAPL-Bucky",
              "JHUAPL-Gecko" = "JHUAPL-Gecko",
              "JHUIDD" = "JHU_IDD-CovidSP",
              "Karlen" = "Karlen-pypm",
              "LANL" = "LANL-GrowthRate",
              "UCSB-ACTS" = "UCSB-ACTS",
              "UCLA" = "UCLA-SuEIR",
              "SI_kJalpha" = "USC-SI_kJalpha"
        )
         )
}

#' Check if estimates vary by location and target_end date
#'
#' Function will return boolean indicator as to whether or not a set of hospitalization data
#' by group (i.e. quantiles functions) varying in number over location and date
#' @param hosp The raw hospitalization data
#' @param locvar string name of location variable (default = "location")
#' @param datevar string name of date variable (default = "target_end_date"); could use "target"
#' @return boolean indicator; if TRUE if the dates are consistent, else false
#' @export
#' @examples
#' estimates_consistent(hosp)
estimates_consistent <- function(hosp, locvar="location", datevar="target_end_date") {
  return(hosp[,.N,by=c(locvar,datevar)][,.("unique_est_per_date" = data.table::uniqueN(N)),by=locvar][, max(unique_est_per_date)]==1)
}

#' Show locations where the number of estimates vary by location and target_end date
#'
#' Function will return a list of locations, along with the number of estimates at each
#' date, if the number of estimates varies over the date range
#' @param hosp The raw hospitalization data
#' @param locvar string name of location variable (default = "location")
#' @param datevar string name of date variable (default = "target_end_date"); could use "target"
#' @param locations vector of character location values (default is NULL); if this is provided, the function
#' will return a list of data tables; one for each location, and each one containing the list
#' of the modeling groups, and the start and end date (and number of days) of each group
#' @return list of locations (if `locations` is NULL) or a list of data tables showing the inconsistencies; one
#' data table per value in `locations`
#' @export
#' @examples
#' show_inconsistent(hosp)
show_inconsistent <- function(hosp, locvar="location", datevar="target_end_date", locations=NULL) {
  if(is.null(locations)) {
    #get the inconsistent locations, and return list of such locations
    inconsistent_locations = hosp[,.N, by=c(locvar,datevar)][,.("est" = data.table::uniqueN(N)),by=locvar][est>1,get(locvar)]
    return(inconsistent_locations)
  } else {
    inconsistent_locations <- hosp[get(locvar) %in% locations,.(start = min(get(datevar)),end = max(get(datevar))),by=c(locvar,"src")]
    if(is.integer(hosp[[datevar]])) {
      inconsistent_locations[,diff:=end-start+1]
    }
    return(split(inconsistent_locations,by=locvar))
  }

}

#' Pull in location information for the location codes in COVID hub processed data
#'
#' This function takes any data frame with a location code and a source for the lookup table. The
#' lookup can either be a path to a csv file, which the function then attempts to read in using
#' `data.table::fread`, or can be data.frame/data.table. Both must have a common variable on
#' which the merge will be based; the name of this common variable can indicated in "location_var"
#' which has default value "location"
#' @param df data frame/data.table for which location information is required
#' @param location_variable character string that indicates the variable in df on which to match
#' the covid19hosp internal location dataset
#' a data table /data.frame that will be the source of the location information
#' @param keepvars an optional character vector of column names to retain after merge
#' @export
#' @return Returns a df merged with the location information
#' @examples
#' add_location_information(ensemble,"../covid19-forecast-hub/data-locations/locations.csv"), location_variable="location")
add_location_info <- function(df,location_variable="location", keepvars=NULL) {

  result <- location_dataset[df,on=location_variable]

  if(!is.null(keepvars)) {
    result <- result[,..keepvars]
  }
  return(result)
}







