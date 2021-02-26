# covid19hosp

### Overview
This is an R package for estimating COVID-19 hospitalization ensemble forecasts from indiviual component models. The package depends on access to a local clone of the [COVID-19 Forecast Hub](https://github.com/reichlab/covid19-forecast-hub)

### Ensembles
The main purpose of the package is to estimate hospitalization ensemble forecasts from the models available at the above repo. There are four possible ensembles that are created. 
- Simulated Linear Pooling Ensemble
- Pre-Smoothed Linear Pooling Ensemble
- Quantile-by-Quntile simple ensemble

### General Approach

- Pull (from the above repository) the most recent contribution from the individual modeling groups, filtering on those groups that include estimates of incident hospitalization.

- Retain only the date/location-specific quantile functions estimates from the modeling groups that provide a full set of 28 daily forward forecasts

- Estimate two simple per-quantile ensembles (vincentization), namely the per-quantile median, and per-quantile mean. 

- Estimate the "Pre-Smoothed Simulated Linear Pooling" ensemble. 
  - This is done by first smoothing the submitted quantile functions for each location, modeling group, and quantile level across the 28 day forecasting period. The smoothing is done using a locally weighted regression smoother with span of 0.75. 
  - In the next step, a penalized cubic spline regression model is used to estimate the continuous quantile function for every location, modeling group, and date.  
  - For each of these regression models, 1000 random draws from the uniform distribution are then applied, to simulate the original model's forecast distribution for that date/location.  
  - The distributions are summed together (i.e. linearly pooled) across the modeling groups to produce a (multi-modal) joint distribution for each date/location, from which we can estimate the ensemble quantiles.


- A "Simulated Linear Pooling" ensemble is similarly estimated, as above, without the pre-smoothing step