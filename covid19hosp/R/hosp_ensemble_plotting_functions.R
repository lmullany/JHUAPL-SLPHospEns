#' Generate a labeled ensemble plot, with individual model estimates
#'
#' Function will receive a data frame with both individual models and the ensemble
#' and return a plot with some labels
#' @param ensemble This is the ensemble data (i.e. produced by generate_slp_ensemble, for example)
#' @param indivmodesl This is the data table containing the hospitalization data (i.e. produced by hosp_processing)
#' @param location_var String variable that indicates the name of the variable that indicates location (default is "location")
#' @param location string variable indicating the location that is being plotting (for plot title purposes only)
#' @param type string variable that describes what type of ensemble is being plotted
#' @param intervals vector of two interval sizes (default 50 and 95)
#' @return returns a plot object
#' @export
#' @examples
#' forecast_plot(ensdata, hospdata,location="US",type="SLP")
#' forecast_plot(ensdata, hospdata)
forecast_plot <- function(ensemble, indivmodels, location_var="location", location=NULL, type=NULL, intervals=c(50,95)) {

  intervals =c(50,95)
  ensemble_input = data.table::copy(ensemble)
  indivmodels_input = data.table::copy(indivmodels)

  src <- data.table::rbindlist(list(
    ensemble_input[,src:="Ensemble"],
    indivmodels_input),use.names = T,fill=TRUE)[,.SD, .SDcols = c("src",location_var, "target_end_date", "quantile", "hospitalizations")]

  if(uniqueN(src[[location_var]])>1) {
    print(unique(src[[location_var]]))
    stop("Location variable in input frames varies; please only submit frames for a single location",call. = F)
  }
  src[[location_var]] <- NULL

  targ_quantiles <- as.character(c(1-(1-intervals/100)/2, (1-intervals/100)/2,0.5))
  targ_quantiles <- targ_quantiles[order(targ_quantiles)]

  src[,quantile:=as.character(quantile)]

  src <- src[quantile %chin% targ_quantiles]
  src <- data.table::dcast(src,src+target_end_date~quantile, value.var="hospitalizations")


  #Get the list of contributing models, for the caption
  modavail <- make_model_caption(src[src!="Ensemble", unique(src)])


  #Get some date values
  mindate <- as.Date(min(src$target_end_date))
  maxdate <- as.Date(max(src$target_end_date))
  date_length <- (as.integer(maxdate-mindate)+1) %/% 7
  dateticks <- seq.Date(mindate, maxdate, length.out=date_length)

  #print(dateticks)
  color_values = c("grey","red")
  color_value_label_ensemble <- ifelse(is.null(type),"Ensemble",paste0("Ensemble (", type, ")"))
  names(color_values) <- c("Individual Model Estimates", color_value_label_ensemble)

  plt <- ggplot2::ggplot(src, ggplot2::aes(target_end_date, `0.5`, group=src)) +
    ggplot2::geom_line(data=src[src!="Ensemble"],ggplot2::aes(color="Individual Model Estimates"), size=1.2) +
    ggplot2::geom_line(data=src[src=="Ensemble"],ggplot2::aes(color=color_value_label_ensemble), size=1.7) +

    ggplot2::geom_ribbon(ggplot2::aes(ymin=`0.025`, ymax=`0.975`, fill="95% Prediction Interval"),
                data=src[src=="Ensemble"],alpha=0.4)+
    ggplot2::geom_ribbon(ggplot2::aes(ymin=`0.25`, ymax=`0.75`, fill="50% Prediction Interval"),
                data=src[src=="Ensemble"],alpha=0.6)+

    ggplot2::xlab("Date") + ggplot2::ylab("Forecasted Hospitalization") +
    ggplot2::ggtitle(paste0("Forecasted Hospitalizations: ", location)) +
    ggplot2::labs(fill="", color="", caption=modavail) +
    ggplot2::scale_x_date(breaks=dateticks) +
    ggplot2::scale_fill_manual(values=c("95% Prediction Interval"="pink", "50% Prediction Interval" = "red")) +
    ggplot2::scale_color_manual(values=color_values) +
    ggplot2::geom_text(data=src[src=="Ensemble"],
              mapping = ggplot2::aes(label=ifelse(target_end_date %in% as.integer(dateticks),
                                         as.character(round(`0.5`,0)),'')),
              hjust=.3,vjust=-.6, size=2, color="blue") +

    ggplot2::geom_text(data=src[src=="Ensemble"],
              mapping = ggplot2::aes(y=`0.025`,label=ifelse(target_end_date %in% as.integer(dateticks),
                                                   as.character(round(`0.025`,0)),'')),
              hjust=.3,vjust=-.6, size=2, color="blue") +
    ggplot2::geom_text(data=src[src=="Ensemble"],
              mapping = ggplot2::aes(y=`0.975`,label=ifelse(target_end_date %in% as.integer(dateticks),
                                                   as.character(round(`0.975`,0)),'')),
              hjust=.3,vjust=-.6, size=2, color="blue") +

    ggplot2::geom_text(data=src[src=="Ensemble"],
              mapping = ggplot2::aes(y=`0.25`,label=ifelse(target_end_date %in% as.integer(dateticks),
                                                  as.character(round(`0.25`,0)),'')),
              hjust=.3,vjust=-.6, size=2, color="blue") +
    ggplot2::geom_text(data=src[src=="Ensemble"],
              mapping = ggplot2::aes(y=`0.75`,label=ifelse(target_end_date %in% as.integer(dateticks),
                                                  as.character(round(`0.75`,0)),'')),
              hjust=.3,vjust=-.6, size=2, color="blue") +

    ggplot2::theme(legend.position="bottom",
          plot.caption = ggplot2::element_text(color = "black", size=8, face="bold")
    ) +
    ggplot2::guides(color=ggplot2::guide_legend(order=1),fill=ggplot2::guide_legend(order=2))

  for(d in dateticks[-1]) {
    plt <- plt + ggplot2::geom_vline(xintercept = d, linetype=2)
  }

  return(plt)
}

#' Generate all plots given an ensemble, and hospitalization frames
#'
#' Function will receive a data frame with both individual models and the ensemble
#' and, for each location in the combined frames, will call the forecast plot function.
#' This is simply a convenience function around the above
#' @param ensemble This is the ensemble data (i.e. produced by generate_slp_ensemble, for example)
#' @param indivmodels This is the data table containing the hospitalization data (i.e. produced by hosp_processing)
#' @param location_var String variable name found in both ensemble and indivmodels that splits by location, and will be used as the title
#' @param type string variable that describes what type of ensemble is being plotted (for labeling only, default is NULL)
#' @param intervals vector of two interval sizes (default 50 and 95)
#' @return returns a frame of plots
#' @export
#' @examples
#' forecast_plot(data,"US","SLP")
plots_by_location <- function(ensemble, indivmodels, location_var, type=NULL, intervals=c(50,95)) {

  locations <- unique(ensemble[[location_var]])
  data.table(location=locations,
             plots = lapply(locations,function(x) {
               forecast_plot(ensemble[get(location_var)==x],indivmodels[get(location_var)==x], location_var = location_var, location = x,type=type, intervals=intervals)
             })
  )
}

#' Generate a plot of individual modeling group raw data
#'
#' This function will return a plot of the modeling group specific median estimates; the inputted
#' data frame should only recieve a single location; model groups can be optionally excluded
#' @param hosp data frame of hospitalization data; preferably filtered down to a single location
#' @param target_var string name of the hospitalization field (default is "value")
#' @param exclude_models optional string vector of names of source models to exclude
#' @param plt_title option plot title
#' @param legendpos optional legend position must be one of "bottom" (default),"right", "left", or "top"
#' @param legrows optional number of legend rows for models (default will be determined by ggplot)
#' @return plot of the individual modeling group's median estimates of the forecasted hospitalizations
#' @export
#' @examples
#' state_plot(hosp[location==24], target_var="hospitalizations")
#' state_plot(hosp[location==24], target_var="value", exclude_models=c("CU-scenario_high"), plt_title("Maryland"))
state_plot <- function(hosp,target_var="value", exclude_models=NULL, plt_title=NULL, legendpos = c("bottom","right","left","top"),legrows=NULL) {
  legendpos <- match.arg(legendpos)

  gginput <- data.table::copy(hosp)[quantile==0.5 & !src %chin% exclude_models][
    ,.SD,.SDcols = c("src","target_end_date","quantile",target_var)]

  setnames(gginput, new=c("src","date","quantile","est"))

  plt <- ggplot2::ggplot(gginput, ggplot2::aes(date, est, color=src, fill=src)) +
    ggplot2::geom_point(size=2) +
    ggplot2::geom_line(size=1.3) +
    ggplot2::scale_color_viridis_d() +
    #geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    ggplot2::theme(legend.position=legendpos) +
    ggplot2::labs(x="",y="Forecasted Hospitalizations", color="Modeling Group", fill="Modeling Group") +
    ggplot2::guides(color=ggplot2::guide_legend(title.position = "top",title.hjust = 0.5,nrow=legrows)) +
    ggplot2::ggtitle(plt_title)

  return(plt)
}

#' Helper function to create contributing models caption, with wrapping
#'
#' Function takes a list of modelnames and returns a caption, with
#' line breaks very `perline` components
#' @param modelname character vector of modelnames
#' @param perline integer value of when to make line breaks, i.e. every `perline`; default
#' is 6
#' @return a character caption
#' @keywords internal
#' @examples
#' make_model_caption(mymodels,4)
make_model_caption <- function(modelnames,perline=6) {
  #get length
  mnl <- length(modelnames)
  #add commas (but not to the last)
  modelnames <- c(paste0(modelnames[1:(mnl-1)], ", "),modelnames[mnl])
  # add title
  modelnames <- c("Contributing Models: ", modelnames)

  newmodelnames <- lapply(seq(mnl+1), function(x) {
    if(x %% perline == 0 & x!=(mnl+1)) return(paste(modelnames[[x]],"\n", collapse=""))
    return(modelnames[[x]])
  })
  return(paste(newmodelnames,collapse=""))
}

