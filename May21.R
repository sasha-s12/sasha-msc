
## Does github work?

set.seed(123)

library(predtools)
library(dplyr)
library(ggplot2)
library(patchwork)

z <- 0.02 #Risk Threshold 

N <- 1000 #Sample size of the future study 

data(gusto)

gusto$kill <- (as.numeric(gusto$Killip)>1)*1 ## Killip Class Binary (0 = KC 1) 
gusto$Y <- gusto$day30 ## Outcome 

## Coding the treatment groups ##
gusto$trt <- ifelse(gusto$tx == "tPA", 1, 0) ## 1 for TPA 0 for other (SKA or SKA combo)
gusto$trt <- factor(gusto$trt)

data_us <- gusto[gusto$regl %in% c(1, 7, 9, 10, 11, 12, 14, 15),] 
data_other <- gusto[!gusto$regl %in% c(1, 7, 9, 10, 11, 12, 14, 15),]

dev_data <- data_other ## training data 
#n <- 500 #Size of the current sample 
#val_data <- data_us[sample(1:(dim(data_us)[1]),n,F),] ## validation data 
val_data <- data_us
val_data$trt_true <- val_data$trt

##### Creating the TBP ####

full_model <- glm(Y ~ trt + age + miloc + pmi + kill + pmin(sysbp,100) + pulse + 
                    trt:(age + miloc + pmi + kill + pmin(sysbp,100) + pulse), data = dev_data, 
                  family=binomial(link="logit"))

drop1(full_model)

## Questions: 
## Why do we include the main effects in the TBP? I understand the reason why we include 
## the interactions is because we want to see how the treatment effects the individual characteristics
## that's what makes it personalized medicine?? 

TBP_model <- glm(Y ~ trt + age + miloc + pmi + kill + pmin(sysbp,100) + pulse + 
                    trt:(pulse), data = dev_data, 
                  family=binomial(link="logit"))

drop1(TBP_model)

#### AIC suggests that the interaction between treatment and pulse is significant 
#### There are probably some other clinically relevant predictors which I will figure out? 

##### Making some PREDICTIONS!!! ####

# Making counterfactual datasets: 

val_data_SKA <- val_data ## coded as 0
val_data_SKA$trt <- factor(0, levels = c(0, 1))

val_data_TPA <- val_data ## coded as 1
val_data_TPA$trt <- factor(1, levels = c(0, 1))

# predicting the 30-day mortality 
pred_mort_SKA <- predict(TBP_model, newdata = val_data_SKA, type = "response") 
## we do type = "response" to get it into probability scale right?

pred_mort_TPA <- predict(TBP_model, newdata = val_data_TPA, type = "response")

# predicting the CATE/TBP 
val_data$TBP <- pred_mort_SKA - pred_mort_TPA ### E[Y(0) - Y(1)]
val_data$pred_mort_SKA <- pred_mort_SKA
val_data$pred_mort_TPA <- pred_mort_TPA

# Creating the new data-frame to hold all these guys: 

TBP_eval <- val_data %>%
  select(
    treatment_received = trt_true, # group (control/treatment)
    observed_event = Y,     # mortality (0 or 1)
    risk_TPA = pred_mort_TPA, # predicted risk if treated
    risk_SKA = pred_mort_SKA, # predicted risk if not treated
    predicted_benefit = TBP  # risk reduction = SKA - TPA
  )

mean(TBP_eval$predicted_benefit)

#### Some clinical way to determine what the benefit threshold would be? ####

#### Lets say D = 0.001
range(TBP_eval$predicted_benefit)
TBP_eval$treat_model <- ifelse(TBP_eval$predicted_benefit > 0.001 , 1, 0)

### calculating the net benefit for congruent/ individuals ####

# Now to find congruent individuals 
### this method seems to be biased??? losing data?? why do we even do this?? 

TBP_eval$congruent <- TBP_eval$treat_model == TBP_eval$treatment_received

# Determine the total number of patients with congruent treatment recommendations

congruent <- TBP_eval %>% filter(congruent == TRUE)

n_congruent <- nrow(congruent)

# the number of these who have an event 
n_congruent_event <- sum(congruent$observed_event)

# the number who are treated
n_congruent_treat <- sum(as.numeric(congruent$treatment_received))


#### nevermind let's just assume average observed event rate among untreated via model####

# Average observed event rate among untreated
untreated_event_pred <- mean(TBP_eval$observed_event[TBP_eval$treat_model == 0] == 1)

treated_event_pred <- mean(TBP_eval$observed_event[TBP_eval$treat_model == 1] == 1)

n0_pred <- sum(TBP_eval$treat_model == 0)
n1_pred <- sum(TBP_eval$treat_model == 1)

NB_pred <- untreated_event_pred - treated_event_pred - (n1_pred*0.001)/(n0_pred+n1_pred)

#### What if we treat all? ####

untreated_event_trtall <- mean(TBP_eval$observed_event[TBP_eval$treatment_received == 0] == 1)

treated_event_trtall <- mean(TBP_eval$observed_event[TBP_eval$treatment_received == 1] == 1)

NB_trtall <- untreated_event_trtall - treated_event_trtall - 0.001

# Whats the advantage? 
NB_trtall - NB_pred

### Function to calculate the net benefit at different thresholds ###

net_benefit <- function(data, TBP_threshold) {
  
  results <- data.frame(
    threshold = numeric(),
    NB_trtall = numeric(),
    NB_pred = numeric(),
    advantage = numeric()
  )
  
  for (i in TBP_threshold) {
    # Model treatment recommendation
    data$treat_model <- ifelse(data$predicted_benefit > i, 1, 0)
    
    # Average observed event rate among untreated (predicted)
    untreated_event_pred <- mean(data$observed_event[data$treat_model == 0])
    
    # Average observed event rate among treated (predicted)
    treated_event_pred <- mean(data$observed_event[data$treat_model == 1])
    
    # Counts
    n0_pred <- sum(data$treat_model == 0)
    n1_pred <- sum(data$treat_model == 1)
    
    # Net benefit for prediction model
    NB_pred <- untreated_event_pred - treated_event_pred - (n1_pred * i) / (n0_pred + n1_pred)
    
    #### Treat all strategy 
    
    # Average observed event rate among untreated (actual)
    untreated_event_trtall <- mean(data$observed_event[data$treatment_received == 0])
    
    # Average observed event rate among treated (actual)
    treated_event_trtall <- mean(data$observed_event[data$treatment_received == 1])
    
    # Net benefit for treat all strategy
    NB_trtall <- untreated_event_trtall - treated_event_trtall - i
    
    # Advantage of prediction model over treat all
    advantage <- NB_pred - NB_trtall
    
    # Add to results
    results <- rbind(results, data.frame(
      threshold = i,
      NB_trtall = NB_trtall,
      NB_pred = NB_pred,
      advantage = advantage
    ))
  }
  
  return(results)
}

TBP_threshold <- seq(0.001, 0.045, by = 0.001)
nb_results <- net_benefit(TBP_eval, TBP_threshold)
print(nb_results)

TBP_threshold_Zoomed <- seq(0.001, 0.00125, by = 0.000001)
nb_results_Zoomed <- net_benefit(TBP_eval, TBP_threshold_Zoomed)
print(nb_results_Zoomed)

p_all <-
  ggplot(nb_results, aes(x = threshold)) +
  geom_line(aes(y = NB_pred, color = "Model Net Benefit")) +
  geom_line(aes(y = NB_trtall, color = "Treat All")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Net Benefit", x = "Threshold", color = "Strategy") 


p_zoomed <- ggplot(nb_results_Zoomed, aes(x = threshold)) +
  geom_line(aes(y = NB_pred, color = "Model Net Benefit")) +
  geom_line(aes(y = NB_trtall, color = "Treat All")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Net Benefit", x = "Threshold", color = "Strategy") 

p_all + p_zoomed






#### Mohsen code (2025.05.21) ####

calc_nb <- function(T)
  
{
  
  TBP_eval$t <- as.numeric(TBP_eval$treatment_received) # RCT assignment 
  
  p_treat <- mean(TBP_eval$t) # proportion of people treated following RCT assginment 
  
  ATE <- sum(TBP_eval$t*TBP_eval$observed_event)/sum(TBP_eval$t) - sum((1-TBP_eval$t)*TBP_eval$observed_event)/sum(1-TBP_eval$t) ## E[Y | A = 1] - E[Y | A = 0]
  
  treated <- TBP_eval[which(TBP_eval$predicted_benefit>=T),] ## returns the rows from the satisfied condition that model says should be treated 
  
  ATT <- sum(treated$t * treated$observed_event)/sum(treated$t) # = E[Y | A = 1, Rz = 1] Among Rz = 1, 
  # this gives the proportion who received treatment in the RCT (A = 1) and actually developed the outcome (Y = 1), 
  # divided by total number treated in that group
  
      - sum((1-treated$t)*treated$observed_event)/sum(1-treated$t) # = E[Y | A = 0, Rz = 1]  
  # Among Rz = 1 this gives the proportion who did not receive treatment (A = 0) and developed the outcome, 
  # divided by total number untreated in that group
  
  c(NB_none=0, NB_model=p_treat*(ATT-T), NB_all=ATE-T) 
  
}

Ts <- (0:100)/1000

res <- matrix(NA, nrow=length(Ts), ncol=3)

for(i in 1:length(Ts))
  
{
  
  res[i,] <- calc_nb(Ts[i])
  
}

plot(Ts, res[,2], type='l', col='red'); lines(Ts, res[,3])


### End: Mohsen code

## Trying the plug-in estimator including the congruent people (2025.05.26)

calc_nb_2 <- function(T) {
  
  TBP_eval$t <- as.numeric(TBP_eval$treatment_received)  # RCT assignment (A)
  
  EY_A0 <- mean(TBP_eval$observed_event[which(TBP_eval$t == 0)])  # E[Y | A = 0] average outcome among those not treated in RCT
  
  TBP_eval$Rz <- as.numeric(TBP_eval$predicted_benefit >= T)  # Rz: model recommendation
  
  EY_ARz <- mean(TBP_eval$observed_event[which(TBP_eval$t == TBP_eval$Rz)])  # E[Y | A = Rz]
  
  p_treat <- mean(TBP_eval$Rz)  # Pr(Rz = 1) proportion of people whom the model recommends treatment
  
  NB_model <- EY_A0 - EY_ARz - T * p_treat # Final output using NB_z formula
  
  # NB_all
  EY_A1 <- mean(TBP_eval$observed_event[which(TBP_eval$t == 1)])  # E[Y | A = 1]
  ATE <- EY_A0 - EY_A1
  NB_all <- ATE - T
  
  c(NB_none = 0, NB_model = NB_model, NB_all = NB_all)
}


Ts <- (0:100)/1000

res <- matrix(NA, nrow=length(Ts), ncol=3)

for(i in 1:length(Ts))
  
{
  
  res[i,] <- calc_nb(Ts[i])
  
}

plot(Ts, res[,2], type='l', col='red'); lines(Ts, res[,3])


plot(Ts, res[,2], type = 'l', col = 'red',
     xlim = c(0, 0.01), 
     ylim = c(-0.1, 0.05))     
lines(Ts, res[,3], col = 'black') 










