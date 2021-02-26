#' Get coverage of predictions
#'
#' Function allows for estimating of coverage of xmin,xmax, given empirical y
#' @param pred_min vector with predicted minimum
#' @param pred_max vector with predicted maximum
#' @param empir_obs vector with observed outcomes
#' @param inclusive logical, default = T, coverage is estimated using inclusive bounds, rather than exclusive
#' @return count (numerator), possible (denominator), and proportion of times (numerator/denominator) that empir_obs falls within pred_min and pred_max, inclusive
#' @export
#' @examples
#' coverage_predictions(lower, upper, observed)
#' coverage_predictions(lower, upper, observed, inclusive=F)

coverage_predictions <- function(pred_min, pred_max, empir_obs, inclusive=T) {
  lenmin = length(pred_min[which(!is.na(pred_min))])
  lenmax = length(pred_max[which(!is.na(pred_max))])
  lenobs = length(empir_obs[which(!is.na(empir_obs))])
  within_bounds = sum(data.table::between(empir_obs,pred_min,pred_max,incbounds=inclusive),na.rm = T)

  return(list(
    "within_bound" = within_bounds,
    "total_obs" = lenmin,
    "coverage" = within_bounds/lenmin

    )
  )
}


#' Get Coverage Stats from a Merged Data Frame
#'
#' This function takes a merged data frame, and calculates coverage stats by any column; the format
#' of the merged data frame is important. It must have both hospital and predicted values. The
#' predicted values must be by quantile, with names like `0.5`, etc.. The function will look for
#' these names, and expects to find them
#' @param dat This is the data frame of merged hospitalization predictions and empirical data
#' @param bycols (default is NULL) a string vector of columns to group by
#' @param qnames (default is NULL) a string vector of the quantile column names; when NULL, the function will look for
#' 23 columns in `dat` starting with 0; if not found will return error
#' @return data table of coverage predictions
#' @export
#' @examples
#' get_coverage_stats(dat=dat)
#' get_coverage_stats(dat=dat,bycols=c("horizon","indicator"))

get_coverage_stats <- function(dat, bycols=NULL,qnames=NULL) {
  if(is.null(qnames)) {
    qnames <- names(dat)[stringr::str_starts(names(dat),"0")]
    if(length(qnames)!=23) stop("Error; the input data frame does not have 23 quantile columns; if less, please input them manually")
  }
  qlength = length(qnames)

  rbindlist(lapply(1:((qlength-1)/2),function(y) {
    dat[,coverage_predictions(pred_min = get(qnames[y]), pred_max=get(qnames[qlength+1-y]), hospitalizations), by=bycols][
      ,`:=`(lower=qnames[y],
            upper=qnames[qlength+1-y],
            diff=as.double(qnames[qlength+1-y])-as.double(qnames[y])
      )]
  }))

}

