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

mean_pred_ben <- mean(TBP_eval$predicted_benefit) # is this where we expect to see the treat all be better than the model or something?


calc_nb <- function(T, TBP_eval)
  
  {
  
  TBP_eval$t <- as.numeric(TBP_eval$treatment_received) # RCT assignment 
  
  p_treat <- mean(TBP_eval$t) # proportion of people treated following RCT assginment 
  
  ATE <- (sum(TBP_eval$t*TBP_eval$observed_event)/sum(TBP_eval$t) - sum((1-TBP_eval$t)*TBP_eval$observed_event)/sum(1-TBP_eval$t)) ## E[Y | A = 1] - E[Y | A = 0]
  
  treated <- TBP_eval[which(TBP_eval$predicted_benefit>=T),] ## returns the rows from the satisfied condition that model says should be treated 
  
  ATT <- (sum(treated$t * treated$observed_event)/sum(treated$t) - # = E[Y | A = 1, Rz = 1] Among Rz = 1, 
  # this gives the proportion who received treatment in the RCT (A = 1) and actually developed the outcome (Y = 1), 
  # divided by total number treated in that group
  
  sum((1-treated$t)*treated$observed_event)/sum(1-treated$t)) # = E[Y | A = 0, Rz = 1]  
  # Among Rz = 1 this gives the proportion who did not receive treatment (A = 0) and developed the outcome, 
  # divided by total number untreated in that group
  
  prop_treated <- nrow(treated) / nrow(TBP_eval) ## proportion of people treated at the threshold
  
  c(NB_none=0, NB_model=p_treat*(ATT-T), NB_all=ATE-T, prop_treated = prop_treated) 
  
}

Ts <- seq(-0.0055, 0.0179, by = 0.0001)  # vector of thresholds 

Ts

res <- matrix(NA, nrow=length(Ts), ncol=4)
colnames(res) <- c("NB_None", "NB_Model", "NB_All", "Prop_treated")

for(i in 1:length(Ts))
  {
  
  res[i,] <- calc_nb(Ts[i],TBP_eval)
  
}

### End: Mohsen code

quantile(TBP_eval$predicted_benefit, probs = c(0.05, 0.95))

hist(TBP_eval$predicted_benefit)

#### Attempting Bootstrapped Confidence Bands 05.27.2024 ####

B <- 1000 # number of iterations

Ts <- seq(-0.0055, 0.0179, by = 0.0001)  # vector of thresholds 

nb_boot_model <- matrix(NA, nrow = B, ncol = length(Ts))
nb_boot_all <- matrix(NA, nrow = B, ncol = length(Ts))

for (b in 1:B) {
  
  val_data_boot <- val_data[sample(nrow(val_data), replace = TRUE), ] # Resampling from the number of rows in val_data
  
  # Computing counterfactual predictions
  val_data_SKA <- val_data_boot
  val_data_SKA$trt <- factor(0, levels = c(0, 1))
  
  val_data_TPA <- val_data_boot
  val_data_TPA$trt <- factor(1, levels = c(0, 1))
  
  pred_mort_SKA <- predict(TBP_model, newdata = val_data_SKA, type = "response")
  pred_mort_TPA <- predict(TBP_model, newdata = val_data_TPA, type = "response")
  
  val_data_boot$TBP <- pred_mort_SKA - pred_mort_TPA
  val_data_boot$pred_mort_SKA <- pred_mort_SKA
  val_data_boot$pred_mort_TPA <- pred_mort_TPA
  
  TBP_eval_boot <- val_data_boot %>%
    select(
      treatment_received = trt_true,
      observed_event = Y,
      #risk_TPA = pred_mort_TPA,
      #risk_SKA = pred_mort_SKA,
      predicted_benefit = TBP
    )
  
  # Computing NB_model and NB_all at each threshold
  for (i in 1:length(Ts)) {
    T_val <- Ts[i]
    nb <- calc_nb(T_val,TBP_eval_boot)
    nb_boot_model[b, i] <- nb["NB_model"]
    nb_boot_all[b, i] <- nb["NB_all"]
  }
}

# 95% confidence bands
lower_band_model <- apply(nb_boot_model, 2, quantile, probs = 0.025, na.rm = TRUE)
upper_band_model <- apply(nb_boot_model, 2, quantile, probs = 0.975, na.rm = TRUE)

lower_band_all <- apply(nb_boot_all, 2, quantile, probs = 0.025)
upper_band_all <- apply(nb_boot_all, 2, quantile, probs = 0.975)

 
plot(Ts, res[, 2], type = "l", col = "red", lwd = 2,
     ylab = "Net Benefit", xlab = "Threshold",
     ylim = c(-0.05, 0.045))


#NB_model
lines(Ts, lower_band_model, col = "red", lty = 2)
lines(Ts, upper_band_model, col = "red", lty = 2)

# NB_all 
lines(Ts, res[, 3], col = "blue", lwd = 2)
lines(Ts, lower_band_all, col = "blue", lty = 2)
lines(Ts, upper_band_all, col = "blue", lty = 2)
 

#hist(apply(nb_boot_all, 1, mean)) ## histogram of the net benefit for treat all strategy after bootstrapping 

#hist(apply(nb_boot_model, 1, mean)) ## histogram of the net benefit for treating by model strategy after bootstrapping 

## overlay proportion 

# Net benefit model
plot(Ts, res[, 2], type = "l", col = "red", lwd = 2,
     
     ylab = "Net Benefit", xlab = "Threshold", 
     ylim = c(-0.012, 0.03))

### 95% confidence bands
lines(Ts, lower_band_model, col = "red", lty = 2)
lines(Ts, upper_band_model, col = "red", lty = 2)

# Net Benefit treat all 
lines(Ts, res[, 3], col = "blue", lwd = 2)

# 95% confidence bands
lines(Ts, lower_band_all, col = "blue", lty = 2)
lines(Ts, upper_band_all, col = "blue", lty = 2)


## Proportion overlay
par(new = TRUE)
plot(Ts, res[, 4], type = "l", col = "darkgreen", lwd = 2,
     axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1))
axis(4)







  
  