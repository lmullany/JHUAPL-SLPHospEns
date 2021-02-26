## ---- include = FALSE-------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=7, 
  fig.height=4)
options(width = 120)


## ----setup, include=FALSE, echo=FALSE---------------------------------------------------------------------------------
# Load necessary libraries for this demonstration
library(psSLP)
library(data.table)
library(ggplot2)

## ----generate_data----------------------------------------------------------------------------------------------------
# Set a sequence of probability levels, p
forecast_1 = data.table::data.table(
  p = c(0.01,0.025,seq(0.05,0.95,0.05),0.975,0.99),
  qp = c(0,1,3,4,6,8,9,10,11,15,17,19,20,25,27,30,31,31,53,82,85,104,164)
)
forecast_1 %>% head()

## ----plot_forecast1, echo=FALSE---------------------------------------------------------------------------------------
plt <- ggplot(data=forecast_1,aes(p,qp)) + 
  geom_point(color="red",size=2) + 
  geom_line(size=1.5) + 
  labs(x="p", y="Q(p)") + 
  ggtitle("Quantile function estimates - Forecast 1")
plt

## ----add_gam_curve, echo=FALSE----------------------------------------------------------------------------------------
plt <- plt + geom_smooth(method="gam",se=F, color='blue', size=1.5)
suppressMessages(print(plt))

## ---------------------------------------------------------------------------------------------------------------------

# Simulate a 1000 observations from the underlying distribution
simulated_1 = generate_slp_predictions(
  x = forecast_1,
  qp_var = "qp",
  p_var = "p")

# plot a histogram of these
ggplot(simulated_1,aes(x=predictions,stat="identity")) + 
  geom_histogram(bins=100,fill="darkgreen") +
  labs(x="Predictions", y="Count") + 
  ggtitle("Simulated distribution from a single quantile function")


## ---------------------------------------------------------------------------------------------------------------------
# add simulate quantiles
forecast_1[, sim_qp:= quantile(simulated_1$predictions, probs = p)]

## ----echo=FALSE-------------------------------------------------------------------------------------------------------
# plot these
ggplot(
  data = melt(forecast_1,
              id.vars = "p",
              measure.vars = c("qp", "sim_qp"),
              value.name = "qp",
              variable.name = "qp_type")[,qp_type:=fifelse(qp_type=="qp", "Original Quantiles", "Estimated from Simulated Distribution")],
  aes(x=p,y=qp,color=qp_type)) + 
  geom_point(size=2)+
  geom_line(size=1.5)+
  labs(x="p", y = "Q(p)", color="Q(p) Type") + 
  theme(legend.position="bottom") + 
  ggtitle("Comparison of estimated and original quantiles from single forecast")


## ---------------------------------------------------------------------------------------------------------------------
# Load the data
data("multiple_models")

## ---------------------------------------------------------------------------------------------------------------------
# Show first 6 rows
multiple_models %>% head()

## ---- echo=F----------------------------------------------------------------------------------------------------------
plt <- ggplot(multiple_models, aes(quantile,outcomes, color=model)) + 
  geom_point() + 
  geom_line() + 
  labs(x="p", y="Q(p)", color="Model") + 
  theme(legend.position="bottom") +
  ggtitle("Multiple quantile functions")

plt

## ---------------------------------------------------------------------------------------------------------------------
ens_vinc <- generate_vinc_ensemble(x=multiple_models,qp_var = "outcomes",p_var = "quantile",pooltype = "mean")

# view the ensemble
ens_vinc %>% head()



## ---- echo=F----------------------------------------------------------------------------------------------------------
plt <- plt + 
  geom_point(data = ens_vinc, aes(x=quantile, y=outcomes, color="Ensemble"), color="black", size=2) + 
  geom_line(data=ens_vinc, aes(x=quantile, y=outcomes, color="Ensemble"), color="black",size=2) + 
  ggtitle("Multiple quantile functions with simple Q-by-Q mean ensemble (black)")

plt


## ---------------------------------------------------------------------------------------------------------------------

# Generate the distributions
sim_distributions <- generate_slp_predictions(multiple_models,qp_var="outcomes",p_var="quantile",byvar="model")

# Now, lets look at the result
sim_distributions[,.N,by=model]

## ---- echo=F----------------------------------------------------------------------------------------------------------
ggplot(sim_distributions, aes(predictions, stat="identity", color=model,fill=model)) + 
  geom_histogram(bins=100)+ 
  facet_wrap(~model, scales="free_y",nrow=2) + 
  theme(legend.position="none") + 
  labs(x="Predictions", y="Count") + 
  ggtitle("Simulated Distributions")
  


## ---- echo=F----------------------------------------------------------------------------------------------------------

vinc_sim_dist <- generate_slp_predictions(ens_vinc,qp_var="outcomes",p_var = "quantile")

ggplot(rbind(vinc_sim_dist[,`:=`(model="vinc",type="Vincentization")], sim_distributions[,type:="SLP Ensemble"]),
       aes(predictions, stat="identity", color=type, fill=type)) +
  geom_density(size=2,alpha=0.3) + 
  theme(legend.position='bottom') +
  labs(x="Predictions", y="Density", color="",fill="") + 
  ggtitle("Joint distribution and distribution simulated from vincentization Q(p)")
  

## ---------------------------------------------------------------------------------------------------------------------
# This line simply estimate the quantiles from the joint (multi-modal) distribution
slp_ensemble <- generate_slp_ensemble(sim_distributions,qp_var="predictions")


## ---- echo=F----------------------------------------------------------------------------------------------------------
gginput <- rbind(slp_ensemble[, type:="SLP Ensemble"],ens_vinc[,type:="Vincentization"], use.names=FALSE)
ggplot(gginput, aes(quantile, predictions, color=type)) + 
  geom_point(size=2) + 
  geom_line(size=1.5) + 
  labs(x="p", y='Q(p)', color="") + 
  ggtitle("Comparison of Q(p) estimates by ensemble approach") + 
  theme(legend.position="bottom")


## ---------------------------------------------------------------------------------------------------------------------
# Load the data
data("multiple_models_over_time")

## ---------------------------------------------------------------------------------------------------------------------

target_models = c("model 1", "model 5", "model 7", "model 8")
target_quantiles = c(0.01, 0.25, 0.5, 0.75, 0.99)

plt <- ggplot(
  data = multiple_models_over_time[
    model %in% target_models
    & quantile %in% target_quantiles], aes(val, outcomes, color=factor(quantile))) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~model, scales="free_y") + 
  theme(legend.position="bottom") + 
  labs(x="Time-Unit", y="Q(p)", color = "p") + 
  ggtitle("Quantile values over time, by model")
plt

## ---------------------------------------------------------------------------------------------------------------------

ps_data <- pre_smooth_quantile_functions(
  x=multiple_models_over_time,
  qp_var="outcomes",
  p_var = "quantile",
  time_var="val",
  byvars="model")


## ---------------------------------------------------------------------------------------------------------------------
ps_data %>% head()

## ---- echo=FALSE------------------------------------------------------------------------------------------------------

target_models = c("model 1", "model 5", "model 7", "model 8")
target_quantiles = c(0.01, 0.25, 0.5, 0.75, 0.99)

plt <- ggplot(
  data = ps_data[model %in% target_models & quantile %in% target_quantiles],
  aes(val, outcomes_hat, color=factor(quantile))) +
  
  geom_point() + 
  geom_line() + 
  facet_wrap(~model, scales="free_y") + 
  theme(legend.position="bottom") + 
  labs(x="Time-Unit", y="Q(p)", color = "p") + 
  ggtitle("Pre-Smoothed quantile values over time, by model")

plt

## ---------------------------------------------------------------------------------------------------------------------
# Get simulate linear pooling of the pre-smoothed quantile functions (by model and time value)
ps_slp <- generate_slp_ensemble(
  generate_slp_predictions(
    x = ps_data, 
    qp_var="outcomes_hat",
    p_var = "quantile",
    byvars = c("model","val")
  ),
  byvars = "val"
)


## ---------------------------------------------------------------------------------------------------------------------

ps_slp_wide <- dcast(
  data = ps_slp[as.character(quantile) %in% c("0.025","0.25","0.5","0.75","0.975")],
  formula = val~quantile,
  value.var = "predictions"
)

ps_slp_wide %>% head()


## ---- echo=FALSE------------------------------------------------------------------------------------------------------
plt <- ggplot(ps_slp_wide, aes(val,`0.5`)) + 
  geom_point() + 
  geom_line() + 
  geom_ribbon(aes(ymin=`0.025`, ymax=`0.975`), alpha=.3) +
  geom_ribbon(aes(ymin=`0.25`, ymax=`0.75`), alpha=.5) + 
  labs(x="Time-unit (val)", y="Q(p)") + 
  ggtitle("Pre-smoothed Simulated Linear Pooling Ensemble")
  
plt

## ---- echo=FALSE, include=FALSE---------------------------------------------------------------------------------------
# Get simulate linear pooling of the raw quantile functions (by model and time value)
slp <- generate_slp_ensemble(
  generate_slp_predictions(
    x = multiple_models_over_time,
    qp_var="outcomes",
    p_var = "quantile",
    byvars = c("model","val")
  ),
  byvars = "val"
)

## ---- echo=FALSE------------------------------------------------------------------------------------------------------
slp_wide <- dcast(
  data = slp[as.character(quantile) %in% c("0.025","0.25","0.5","0.75","0.975")],
  formula = val~quantile,
  value.var = "predictions"
)

plt <- ggplot(slp_wide, aes(val,`0.5`)) + 
  geom_point() + 
  geom_line() + 
  geom_ribbon(aes(ymin=`0.025`, ymax=`0.975`), alpha=.3) +
  geom_ribbon(aes(ymin=`0.25`, ymax=`0.75`), alpha=.5) + 
  labs(x="Time-unit (val)", y="Q(p)") + 
  ggtitle("Simulated Linear Pooling Ensemble (No Smoothing)")
  
plt



