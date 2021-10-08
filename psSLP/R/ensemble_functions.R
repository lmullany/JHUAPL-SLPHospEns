
#' Generate simulated linear pooling predictions
#'
#' This function simulates a set of simulated outputs, given quantile function. Provide a data
#' frame `x` that has at least two columns, the quantile probability (p) and the threshold value
#' Q(p); indicating these as the `p_var` and the `qp_var`.  If you have more than one
#' quantile function set, then indicate the variables that identify these sets in
#' `byvars`)
#'
#' @param x This is the data frame/data.table that contains the `p_var` and `qp_var` and any
#' grouping variables `byvars` (i.e. defining set of quantile functions).
#' @param qp_var This is the name of the variable that holds the estimate of the number
#' of hospitalizations.
#' @param p_var This is the name of the variable that holds the probability, p
#' @param byvars (default = NULL), a character vector indicating
#' the grouping variables. Each unique group should represent a single quantile function (Q(p))
#' @param draw_size (default=1000), an integer value of the number of simulated forecasts to
#' make for each quantile function set
#' @param parallel: string value either "auto","on","off". By default, the function will run in parallel
#' using half the available threads if the number of gam models to estimates exceeds 1000 (this is "auto")
#' mode. Setting this to "on" ("off") will force parallel processing on ("off") regardless of gam models
#' @param threads default NULL, set integer number of threads to use for parallel processing
#' @param verbose (default = TRUE) a boolean indicator to show verbose status updates during calculation
#' @param gam_model_configuration; provide a list of additional parameters to the gam model. Currently, the
#' default list is `list(family="gaussian", k=10)` to control the model distribution and the number of knots
#' @return Returns a data table of predictions (size=`draw_size`) for each grouping variable defined via `byvars`
#' @import data.table
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#' @export
#' @examples
#' generate_slp_predictions(hosp=my_hosp_data, target_variable="value",draw_size=500, byvars=c("location,"source","date))
generate_slp_predictions <- function(x,
                                  qp_var,
                                  p_var,
                                  byvars = NULL,
                                  draw_size = 1000,
                                  parallel = c("auto","on","off"),
                                  threads = NULL,
                                  verbose=T,
                                  gam_model_configuration=list(
                                    family="gaussian",
                                    k=10)) {

  parallel = match.arg(parallel)

  # check that if a gam_model_configuration list has been provided, then
  # any elements that are null get filled with default.

  if(is.null(gam_model_configuration[["family"]])) gam_model_configuration[["family"]] <- "gaussian"
  if(is.null(gam_model_configuration[["k"]])) gam_model_configuration[["k"]] <- 10

  # Check Inputs
  check_inputs(input_list = list("qp_var" = qp_var,"p_var" = p_var,"byvars" = byvars),x = x)

  # A function to generate predictions for each model
  sim_model_output <- function(gam_model, numpredictions) {
    phat <- predict(gam_model,newdata=data.frame(quantile=runif(n=numpredictions)))
    return(as.vector(phat))
  }

  input_df = data.table::setDT(data.table::copy(x))
  input_df[,"quantile":=get(p_var)]

  #if byvars is null, this should be very quick (one model only)
  if(is.null(byvars)) {
    gam_model <- suppressWarnings(
      mgcv::gam(get(qp_var)~s(quantile,bs="cs",
                              k=gam_model_configuration$k),data=input_df, model=F, family=gam_model_configuration$family)
    )
    predictions = sim_model_output(gam_model,numpredictions = draw_size)
    if(gam_model_configuration$family == "poisson") {
      print("poisson_predictions!")
      predictions = exp(predictions)
    }

    return(data.table::data.table("predictions"=predictions))
  }

  nmods <- data.table::uniqueN(input_df,by=byvars)
  parallel_config = get_parallel_status(nmods,parallel,threads)


  if(verbose) {
    if(parallel_config[["parallel"]]){
      cat("Running",nmods,"gam models in parallel.. will take some time\n")
    } else {
      cat("Running",nmods,"gam models.. will take some time\n")
    }

    if(!parallel_config[["parallel"]] & nmods>1000) cat("Not running in parallel mode, but more than 1000 models to be run.. consider setting parallel to auto/on\n")
  }

  #Now (possibly) parallelize the gam model and prediction process, for the (possibly thousands of) groups
  elapsed <- system.time({
    if(parallel_config[["parallel"]]) {
      clusters <- parallel::makeCluster(parallel_config[["threads"]])
      doParallel::registerDoParallel(clusters)
      `%par_type%` <- `%dopar%`
    } else {
      `%par_type%` <- `%do%`
    }
    result <- foreach::foreach(i=split(input_df,by=byvars),.combine = rbind, .verbose=F,.packages="data.table",.inorder = F) %par_type% {
        gam_model <- suppressWarnings(
          mgcv::gam(get(qp_var)~s(quantile,bs="cs",k=gam_model_configuration$k),data=i, model=F, family=gam_model_configuration$family)
        )
        predictions = sim_model_output(gam_model,draw_size)
        if(gam_model_configuration$family == "poisson") {
          predictions = exp(predictions)
        }
        keyv <- lapply(byvars,function(x) unique(i[[x]]))
        res = data.table::data.table(predictions=predictions)
        for(kv in seq(length(keyv))) data.table::set(x=res,j = byvars[kv],value = keyv[kv])
        res[]
    }
    if(parallel_config[["parallel"]]) {
      parallel::stopCluster(clusters)
    }
  })

  if(verbose) cat(paste0("Complete in ",round(elapsed[3],2)," seconds."),"\n")

  return(result)
}

#' Generate simulated linear pooling ensemble
#'
#' This function takes a set of simulated predictions (for example, returned by `generate_slp_predictions`,
#' and returns the estimates of the quantile function Q(p) for the set. Optionally, multiple sets of
#' simulated predictions can be given, with sets distinguished by `byvars` parameter
#'
#' @param slp_predictions This is the set of simulated predictions returned by `generate_slp_predictions`; it
#' must contain a column named "predictions"
#' @param qp_var string name of the resulting combined variable (defaults to "predictions")
#' @param byvars (default = NULL), a character vector indicating the variables that uniquely define
#' a set of predictions
#' @param quantiles a numeric vectors of values between 0 and 1 exclusive; these are the quantile levels to
#' estimate. The default set contains 23 levels (`c(0.01,0.025,seq(0.05,0.95,0.05),0.975,0.99)`)
#' @return Returns a data table of the simulated linear pooling ensemble
#' @export
#' @examples
#' generate_slp_ensemble(slp_predictions = slp)
#' generate_slp_ensemble(slp_predictions, byvars=c("location","date"),quantiles=c(0.025,0.25,0.5,0.75,0.975))

generate_slp_ensemble <- function(slp_predictions,
                                  qp_var = "predictions",
                                  byvars = NULL,
                                  quantiles = c(0.01,0.025,seq(0.05,0.95,0.05),0.975,0.99)) {

  # Check Inputs
  check_inputs(input_list = list("predictions" = "predictions","byvars" = byvars),x = slp_predictions)


  ensemble <- slp_predictions[,.("value" = quantile(predictions,probs = quantiles)), by=byvars][,quantile:=quantiles, by=byvars]
  setnames(ensemble, old="value",new=qp_var)

  return(ensemble[,.SD, .SDcols=c(byvars,"quantile",qp_var)])
}


#' Generate vincentization / quantile-over-quantile ensemble
#'
#' This function takes a set of quantile functions, and returns the mean (default) or median
#' quantile level over the quantiles. This is also known as the vincentization approach.
#'
#' @param x This is the data frame/data.table that contains the `p_var` and `qp_var` and any
#' grouping variables `byvars` (i.e. defining set of quantile functions).
#' @param qp_var This is the name of the variable that holds the estimate of the number
#' of hospitalizations.
#' @param p_var This is the name of the variable that holds the probability, p, `qp_var` will be pooled
#' over equal levels of `p_var`
#' @param byvars (default = NULL), a character vector indicating which variables distinguish different
#' ensembles
#' @param pooltype a string indicating whether "mean" (default) or "median" of the quantile levels should
#' be returned
#' @return Returns a data table of the quantile-by-quantile or vincentization ensemble
#' @export
#' @examples
#' generate_vinc_ensemble(hosp, "value")
#' generate_vinc_ensemble(hosp, target_variable="val", pooltype="median")
#' generate_vinc_ensemble(hosp, target_variable="val", pooltype="median", byvars=c("loc","date","quantile")
generate_vinc_ensemble <- function(x,
                                   qp_var,
                                   p_var,
                                   byvars = NULL,
                                   pooltype=c("mean","median")) {

  # Check Inputs
  check_inputs(input_list = list("qp_var" = qp_var,"p_var" = p_var,"byvars" = byvars),x = x)

  # the pooling type must be median or mean, and use match.fun to get the right function
  pooltype = match.fun(match.arg(pooltype))

  ensemble <- x[,.("value" = pooltype(get(qp_var),na.rm=T)), by=c(byvars,p_var)]
  setnames(ensemble,old="value",new=qp_var)

  return(ensemble[,.SD, .SDcols=c(byvars,p_var,qp_var)])

}

#' Pre-smooth individual model quantile function, prior to applying the SLP method
#'
#' This function takes in a set of individual models, and pre-smooths the quantile
#' function to even out day-to-day volatility, prior to inserting into the SLP
#' prediction process. The function using a locally weight regression smoothing approach
#' @param x the source data table containing the quantile function(s) over some unit, generally time
#' @param qp_var The name of the variable holding the $Q(p)$ values
#' @param p_var The name of the variable in `x` holding the p values
#' @param time_var The name of the variable in `x` holding the unit (i.e. time) over which you want to pre-smooth Q(p) for
#' fixed values of `p_var`
#' @param byvars The vector of additional (default =NULL) variables in `x` that define additional groups
#' @param loess_span numeric value between 0 and 1 (default=0.75) indicating the span of data which will be used in the
#' locally weighted regression process
#' @param smooth_col_name string name of new column that will contain the smoothed Q(p) values; the default is to
#' produce a column name `<qp_var>_hat` (unless it already exists
#' @return Returns the a data table wih `qp_var`,`p_var`,`time_var`,`byvars` (if any), and the new `smooth_col_name`
#' @export
#' @examples
#' pre_smooth_quantile_functions(x=hosp,qp_var="hospitalizations",p_var="quantile", time_var="target_end_date", byvars=c("src","location"), smooth_col_name="value_sm")


pre_smooth_quantile_functions <- function(x, qp_var, p_var, time_var, byvars=NULL, smooth_col_name=paste0(qp_var,"_hat"),loess_span=0.75 ) {

  byvars = c(p_var,byvars)

  # Get the number of unique loess models

  if(!is.null(byvars)) n_loess_models <- uniqueN(x[,..byvars])
  else n_loess_models <- 1

  cat("\nSmoothing",n_loess_models,"quantile levels\n")

  pre_smooth <- x[,list(list(stats::loess(get(qp_var)~as.numeric(get(time_var)), data=.SD,span = loess_span))), by=byvars]

  #add those predicted values back in, but be careful to properly sort
  data.table::setorderv(x,cols = c(byvars, time_var))
  data.table::setorderv(pre_smooth,cols=byvars)



  return(
    data.table::copy(x)[,eval(smooth_col_name):=unlist(lapply(pre_smooth$V1,predict))][,.SD,.SDcols = c(qp_var,time_var,byvars,smooth_col_name)]
  )

}


#' Helper function to check inputs
#'
#' Function for checking varnames passed as strings
#' @param input_list list of inputs variables
#' @param x input dataframe
#' @keywords internal

check_inputs <- function(input_list,x) {

  for(i in names(input_list)) {
    for(el in input_list[[i]]) {
      if(!el %in% names(x)) {
        stop(paste0(i," element `",el,"` not found in input dataframe"),call. = F)
      }
    }
  }
}

#' Helper function to set up parallelization config
#'
#' Function sets up parallelization configuration list
#' @param nmods integer number of models
#' @param parallel string, one of "auto","on","off"
#' @param clusters integer number of threads (or NULL)
#' @keywords internal
get_parallel_status <- function(nmods,parallel,clusters) {

  if(is.character(clusters)) {
    stop("Number of threads must be NULL or numeric",call. = F)
  }

  parallel_config=list(parallel = FALSE,threads=NA)

  #if parallel off, just return the base configuration
  if(parallel=="off" | (parallel=="auto" & nmods<1000)) {
    return(parallel_config)
  }
  parallel_config[["parallel"]] <- TRUE
  # get the minimum of number of cores/clusters
  max_threads = min(parallel::detectCores(),clusters,na.rm = T)
  #Now, set the number of clusters to detect_cores/2
  if(is.null(clusters) || is.na(clusters)) parallel_config[["threads"]] <- max_threads/2
  else parallel_config[["threads"]] <- min(max_threads,clusters,na.rm=T)
  return(parallel_config)
}


