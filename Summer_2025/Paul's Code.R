#Reproducibility 
set.seed(123) 

#Packages used 
library(predtools)
library(dplyr)
library(ggplot2)

##### Setting up the Data #####

data(gusto)

gusto$kill <- (as.numeric(gusto$Killip)>1)*1 ## Killip Class Binary (0 = KC 1) 
gusto$Y <- gusto$day30 ## Outcome 

## Coding the treatment groups ##
gusto$trt <- ifelse(gusto$tx == "tPA", 1, 0) ## 1 for TPA 0 for other (SKA or SKA combo)
gusto$trt <- factor(gusto$trt) 

data_us <- gusto[gusto$regl %in% c(1, 7, 9, 10, 11, 12, 14, 15),] 
data_other <- gusto[!gusto$regl %in% c(1, 7, 9, 10, 11, 12, 14, 15),]

dev_data <- data_other ## training data 
val_data <- data_us ## validation data 
val_data$trt_true <- val_data$trt ## vector trt_true contains the true treatment assignment assigned via RCT

##### Creating the TBP ####

full_model <- glm(Y ~ trt + age + miloc + pmi + kill + pmin(sysbp,100) + pulse + 
                    trt:(age + miloc + pmi + kill + pmin(sysbp,100) + pulse), data = dev_data, 
                  family=binomial(link="logit"))

TBP_model <- glm(Y ~ trt + age + miloc + pmi + kill + pmin(sysbp,100) + pulse + 
                   trt:(pulse), data = dev_data, 
                 family=binomial(link="logit"))

drop1(TBP_model)

##### Making some PREDICTIONS!!! ####

# First we need to figure out what the CATE is, and to do so we need the counterfactual data sets

# Making counter factual data sets: 

# The treatment is not applied 
val_data_SKA <- val_data ## coded as 0
val_data_SKA$trt <- factor(0, levels = c(0, 1))

# The treatment is applied 
val_data_TPA <- val_data ## coded as 1
val_data_TPA$trt <- factor(1, levels = c(0, 1))

# predicting the 30-day mortality:

# treatment is not applied 
pred_mort_SKA <- predict(TBP_model, newdata = val_data_SKA, type = "response") 

# treatment is applied 
pred_mort_TPA <- predict(TBP_model, newdata = val_data_TPA, type = "response")

# predicting the CATE/TBP 
val_data$TBP <- pred_mort_SKA - pred_mort_TPA ### E[Y(0)] - E[Y(1)]
val_data$pred_mort_SKA <- pred_mort_SKA
val_data$pred_mort_TPA <- pred_mort_TPA

# Creating the new data-frame to hold all these guys: 

TBP_eval <- val_data %>%
  select(
    treatment_received = trt_true, # group (control/treatment) This is our "A"
    observed_event = Y,     # mortality (0 or 1) This is our outcome 
    #risk_TPA = pred_mort_TPA, # predicted risk if treated
    #risk_SKA = pred_mort_SKA, # predicted risk if not treated
    predicted_benefit = TBP  # risk reduction = SKA - TPA
  )


summary(TBP_eval)

NB.PG <- function(y, a, pb) {
  ### doing quickly, but ...
  ### should check all three inputs are same length
  ### should check a, y numeric, binary
  ### should check pb is numeric, in [-1,1]
  
  ## outcome prevalence in trial arms
  y.bar.cntrl <- mean(y[a==0])
  
  y.bar.trt <- mean(y[a==1])
  
  ## very fine grid
  t.grd <- sort(unique(TBP_eval$predicted_benefit))
  
  nb.model <- rep(NA, length(t.grd)) ### will be output
  
  for (i in 1:length(t.grd)) {
    ### who is compliant at this threshold?
    
    ndx <- as.logical(a) == (pb > t.grd[i]) ### returns a list of the congruent individuals 
    
    nb.model[i] <- y.bar.cntrl - mean(y[ndx]) - t.grd[i]*mean(pb > t.grd[i])
    
  }
  
  list(t.grd=t.grd, nb.model=nb.model,
       nb.all=y.bar.cntrl-y.bar.trt-t.grd)
}

tmp <- NB.PG(y=TBP_eval$observed_event,
            
             a = as.numeric(as.character(TBP_eval$treatment_received)),
            
            pb=TBP_eval$predicted_benefit)


## nb.model
plot(tmp$t.grd, tmp$nb.model,xlab="threshold",ylab="NB",
     xlim=quantile(TBP_eval$predicted_benefit, c(.01,.99)),
     ylim=c(-0.02, 0.05), type="l")

## nb.all
points(tmp$t.grd, tmp$nb.all,col="blue",type="l")

## nb.nobody
abline(h=0,col="red")














