### Let's try this out on real data!!! 

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

TBP_eval <- data.frame(
  y  = TBP_eval$observed_event,
  a  = as.numeric(as.character(TBP_eval$treatment_received)),
  pb = TBP_eval$predicted_benefit
)
















##### Bringing in the Estimators #####

#### Congruent Based Estimator (Paul's Estimator) ####
NB.CB <- function(y, a, pb, t.grd=NULL) {
  
  ### doing quickly, but 
  ### should check all three inputs are same length 
  ### should check a, y numeric, binary  
  ### should check pb is numeric, in [-1,1]
  
  ## outcome prevalence in trial arms
  y.bar.cntrl <- mean(y[a==0]);  
  y.bar.trt <- mean(y[a==1])
  
  ##  very fine grid
  if(is.null(t.grd))
  {
    t.grd <- sort(unique(pb))
  }
  
  nb.model <- rep(NA, length(t.grd))  ### will be output
  
  for (i in 1:length(t.grd)) {
    
    r <- as.integer(pb > t.grd[i])
    
    E_Y_use_rule <- mean(y[a == 0 & r == 0]) * mean(r == 0) + mean(y[a == 1 & r == 1]) * mean(r == 1)
    
    nb.model[i] <- y.bar.cntrl - E_Y_use_rule - t.grd[i]*mean(pb > t.grd[i])   
    
  }
  
  list(t.grd=t.grd, nb.model=nb.model, 
       nb.all=y.bar.cntrl-y.bar.trt-t.grd)  
} 

### Calculating the Standard Errors for the Congruent Based Estimator ####
compute_se_cb <- function(y, a, pb, t.grd) {
  
  n <- length(y)
  se.pg <- numeric(length(t.grd))
  theta_hat <- numeric(length(t.grd))
  ci.lower <- numeric(length(t.grd))
  ci.upper <- numeric(length(t.grd))
  
  for (i in seq_along(t.grd)) {
    
    r <- as.integer(pb > t.grd[i])
    
    counts <- table(factor(a, levels = 0:1),
                    factor(r, levels = 0:1),
                    factor(y, levels = 0:1))
    
    q1 <- counts["0", "0", "0"] / n
    q2 <- counts["0", "0", "1"] / n
    q3 <- counts["0", "1", "0"] / n
    q4 <- counts["0", "1", "1"] / n
    q5 <- counts["1", "0", "0"] / n
    q6 <- counts["1", "0", "1"] / n
    q7 <- counts["1", "1", "0"] / n
    q8 <- counts["1", "1", "1"] / n
    
    q <- c(q1, q2, q3, q4, q5, q6, q7, q8)
    
    d1 <- q1 + q2 + q3 + q4
    benefit <- (q2 + q4) / d1
    
    harm_control <- (q2 / (q1 + q2)) * (q1 + q2 + q5 + q6)
    harm_treated <- (q8 / (q7 + q8)) * (q3 + q4 + q7 + q8)
    
    p <- q3 + q4 + q7 + q8
    
    theta_hat[i] <- benefit - (harm_control + harm_treated) - t.grd[i] * p
    
    grad <- numeric(8)
    grad[1] <- (q2 / (q1 + q2)) - (q2 * (q1 + q2 + q5 + q6)) / (q1 + q2)^2 - (q2 + q4) / d1^2
    grad[2] <- -((q1^2 + 2 * q1 * q2 + q2^2 + q1 * (q5 + q6)) / (q1 + q2)^2) + 1 / d1
    grad[3] <- -q8 / (q7 + q8) - t.grd[i] - (q2 + q4) / d1^2
    grad[4] <- -q8 / (q7 + q8) - t.grd[i] - (q2 + q4) / d1^2 + 1 / d1
    grad[5] <- -q2 / (q1 + q2)
    grad[6] <- -q2 / (q1 + q2)
    grad[7] <- (q3 * q8 + q4 * q8) / (q7 + q8)^2 - t.grd[i]
    grad[8] <- (q8 * (q3 + q4 + q7 + q8) - t.grd[i] * (q7 + q8)^2 - 
                  (q7 + q8) * (q3 + q4 + q7 + 2 * q8)) / (q7 + q8)^2
    
    Sigma <- matrix(0, 8, 8)
    for (j in 1:8) {
      for (k in 1:8) {
        if (j == k) {
          Sigma[j, k] <- q[j] * (1 - q[j]) / n
        } else {
          Sigma[j, k] <- -q[j] * q[k] / n
        }
      }
    }
    
    var_theta <- t(grad) %*% Sigma %*% grad
    se.pg[i] <- sqrt(as.numeric(var_theta))
    
    ci.lower[i] <- theta_hat[i] - 1.96 * se.pg[i]
    ci.upper[i] <- theta_hat[i] + 1.96 * se.pg[i]
  }
  
  return(list(
    t.grd = t.grd,
    theta_hat = theta_hat,
    se = se.pg,
    ci.lower = ci.lower,
    ci.upper = ci.upper
  ))
}

#### Calculating the Shared and Unshared Bayesian Standard Errors for the Congruent Based Estimator #### 

bayesian_NB_CB_Unshared <- function(y, a, pb, t.grd, m = 1000) {
  
  n <- length(y)
  nb.posterior <- matrix(NA, nrow = m, ncol = length(t.grd))
  nb.cb <- numeric(length(t.grd))     # plug-in
  nb.mean <- numeric(length(t.grd))   # posterior mean
  nb.lower <- numeric(length(t.grd))  # 2.5%
  nb.upper <- numeric(length(t.grd))  # 97.5%
  
  
  for (i in seq_along(t.grd)) {
    
    r <- as.integer(pb > t.grd[i])
    
    counts <- table(factor(a, levels = 0:1),
                    factor(r, levels = 0:1),
                    factor(y, levels = 0:1))
    
    q1 <- counts["0", "0", "0"]
    q2 <- counts["0", "0", "1"]
    q3 <- counts["0", "1", "0"]
    q4 <- counts["0", "1", "1"]
    q5 <- counts["1", "0", "0"]
    q6 <- counts["1", "0", "1"]
    q7 <- counts["1", "1", "0"]
    q8 <- counts["1", "1", "1"]
    
    q_counts <- c(q1, q2, q3, q4, q5, q6, q7, q8)
    alpha <- q_counts + 1
    qhat <- q_counts / sum(q_counts)
    
    # Plug-in estimator
    harm_ctrl <- if ((qhat[1] + qhat[2]) > 0) {
      (qhat[2] / (qhat[1] + qhat[2])) * (qhat[1] + qhat[2] + qhat[5] + qhat[6])
    } else 0
    
    harm_trt <- if ((qhat[7] + qhat[8]) > 0) {
      (qhat[8] / (qhat[7] + qhat[8])) * (qhat[3] + qhat[4] + qhat[7] + qhat[8])
    } else 0
    
    benefit <- (qhat[2] + qhat[4]) / (qhat[1] + qhat[2] + qhat[3] + qhat[4])
    p <- qhat[3] + qhat[4] + qhat[7] + qhat[8]
    
    nb.cb[i] <- benefit - (harm_ctrl + harm_trt) - t.grd[i] * p
    
    # Posterior draws from Dirichlet
    gamma_draws <- matrix(
      rgamma(m * 8, shape = rep(alpha, each = m)),
      nrow = m, ncol = 8, byrow = FALSE
    )
    dirichlet_draws <- gamma_draws / rowSums(gamma_draws)
    
    for (j in 1:m) {
      q <- dirichlet_draws[j, ]
      
      harm_ctrl_j <- if ((q[1] + q[2]) > 0) {
        (q[2] / (q[1] + q[2])) * (q[1] + q[2] + q[5] + q[6])
      } else 0
      
      harm_trt_j <- if ((q[7] + q[8]) > 0) {
        (q[8] / (q[7] + q[8])) * (q[3] + q[4] + q[7] + q[8])
      } else 0
      
      benefit_j <- (q[2] + q[4]) / (q[1] + q[2] + q[3] + q[4])
      p_j <- q[3] + q[4] + q[7] + q[8]
      
      nb.posterior[j, i] <- benefit_j - (harm_ctrl_j + harm_trt_j) - t.grd[i] * p_j
    }
    
    nb.mean[i] <- mean(nb.posterior[, i], na.rm = TRUE)
    nb.lower[i] <- quantile(nb.posterior[, i], 0.025, na.rm = TRUE)
    nb.upper[i] <- quantile(nb.posterior[, i], 0.975, na.rm = TRUE)
  }
  
  list(
    t.grd = t.grd,
    nb.plugin = nb.cb,
    nb.posterior = nb.posterior,
    nb.posterior.mean = nb.mean,
    ci.lower = nb.lower,
    ci.upper = nb.upper
  )
}

##### Calculating the bootstrap standard errors for the congruent based estimator ####

bootstrap_NB_CB <- function(y, a, pb, t.grd = NULL, B = 1000) {
  
  n <- length(y)
  
  if (is.null(t.grd)) {
    t.grd <- sort(unique(pb))
  }
  
  nb.boot <- matrix(NA, nrow = B, ncol = length(t.grd))
  
  for (b in 1:B) {
    
    idx <- sample(1:n, size = n, replace = TRUE)
    
    boot.out <- NB.CB(y = y[idx],
                      a = a[idx],
                      pb = pb[idx],
                      t.grd = t.grd)
    
    nb.boot[b, ] <- boot.out$nb.model
  }
  
  nb.lower <- apply(nb.boot, 2, quantile, probs = 0.025, na.rm = TRUE)
  nb.upper <- apply(nb.boot, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  list(
    t.grd = t.grd,
    nb.boot = nb.boot,
    
    ci.lower = nb.lower,
    ci.upper = nb.upper
  )
}

#### ATT Based Estimator ####
NB.ATTB <- function(y, a, pb, t.grd=NULL) {
  
  ### doing quickly, but ...
  ### should check all three inputs are same length 
  ### should check a, y numeric, binary  
  ### should check pb is numeric, in [-1,1]
  
  ## outcome prevalence in trial arms
  y.bar.cntrl <- mean(y[a==0]);  y.bar.trt <- mean(y[a==1])
  
  ##  very fine grid
  if(is.null(t.grd))
  {
    t.grd <- sort(unique(pb))
  }
  
  nb.model <- rep(NA, length(t.grd))  ### will be output
  
  for (i in 1:length(t.grd)) {
    
    ### who is above threshold?
    ndx <- which(pb >= t.grd[i])
    nb.model[i] <- (mean(y[ndx][a[ndx]==0])-mean(y[ndx][a[ndx]==1])-t.grd[i])*mean(pb > t.grd[i]) 
  }
  
  list(t.grd=t.grd, nb.model=nb.model, 
       nb.all=y.bar.cntrl-y.bar.trt-t.grd)  
}  

compute_se_attb <- function(y, a, pb, t.grd) {
  
  n <- length(y)
  se.ms <- numeric(length(t.grd))
  theta_hat <- numeric(length(t.grd))
  ci.lower <- numeric(length(t.grd))
  ci.upper <- numeric(length(t.grd))
  
  for (i in seq_along(t.grd)) {
    
    r <- as.integer(pb >= t.grd[i])
    
    counts <- table(factor(a, levels = 0:1),
                    factor(r, levels = 0:1),
                    factor(y, levels = 0:1))
    
    q1 <- counts["0", "0", "0"]
    q2 <- counts["0", "0", "1"]
    q3 <- counts["0", "1", "0"]
    q4 <- counts["0", "1", "1"]
    q5 <- counts["1", "0", "0"]
    q6 <- counts["1", "0", "1"]
    q7 <- counts["1", "1", "0"]
    q8 <- counts["1", "1", "1"]
    
    q <- c(q1, q2, q3, q4, q5, q6, q7, q8) / n
    
    p <- q[3] + q[4] + q[7] + q[8]
    EY_r1_a0 <- q[4] / (q[3] + q[4])
    EY_r1_a1 <- q[8] / (q[7] + q[8])
    theta_hat[i] <- (EY_r1_a0 - EY_r1_a1 - t.grd[i]) * p
    
    grad <- numeric(8)
    grad[1] <- 0
    grad[2] <- 0
    grad[3] <- (-q[4] / (q[3] + q[4])^2) * p + (EY_r1_a0 - EY_r1_a1 - t.grd[i])
    grad[4] <- (q[3] / (q[3] + q[4])^2) * p + (EY_r1_a0 - EY_r1_a1 - t.grd[i])
    grad[5] <- 0
    grad[6] <- 0
    grad[7] <- (q[8] / (q[7] + q[8])^2) * p + (EY_r1_a0 - EY_r1_a1 - t.grd[i])
    grad[8] <- (-q[7] / (q[7] + q[8])^2) * p + (EY_r1_a0 - EY_r1_a1 - t.grd[i])
    
    Sigma <- matrix(0, 8, 8)
    for (j in 1:8) {
      for (k in 1:8) {
        if (j == k) {
          Sigma[j, k] <- q[j] * (1 - q[j]) / n
        } else {
          Sigma[j, k] <- -q[j] * q[k] / n
        }
      }
    }
    
    var_theta <- t(grad) %*% Sigma %*% grad
    se.ms[i] <- sqrt(as.numeric(var_theta))
    
    # Confidence intervals
    ci.lower[i] <- theta_hat[i] - 1.96 * se.ms[i]
    ci.upper[i] <- theta_hat[i] + 1.96 * se.ms[i]
  }
  
  return(list(
    t.grd = t.grd,
    theta_hat = theta_hat,
    se = se.ms,
    ci.lower = ci.lower,
    ci.upper = ci.upper
  ))
}

#### Calculating the Bayesian Standard Errors for the ATT-Based Estimator #### 

#### Calculating the Bayesian Standard Errors for the ATT-Based Estimator Unshared ####
bayesian_NB_ATTB_Unshared <- function(y, a, pb, t.grd, m = 1000) {
  
  nb.posterior <- matrix(NA, nrow = m, ncol = length(t.grd))
  nb.ms <- numeric(length(t.grd))          # Plug-in estimates
  nb.mean <- numeric(length(t.grd))        # Posterior mean for each t
  nb.low <- numeric(length(t.grd))         # Lower 2.5% percentile
  nb.high <- numeric(length(t.grd))        # Upper 97.5% percentile
  
  for (i in seq_along(t.grd)) {
    
    r <- as.integer(pb >= t.grd[i])
    
    # Contingency table (A, R, Y)
    counts <- table(factor(a, levels = 0:1),
                    factor(r, levels = 0:1),
                    factor(y, levels = 0:1))
    
    # Named extraction for q1 to q8
    q1 <- counts["0", "0", "0"]
    q2 <- counts["0", "0", "1"]
    q3 <- counts["0", "1", "0"]
    q4 <- counts["0", "1", "1"]
    q5 <- counts["1", "0", "0"]
    q6 <- counts["1", "0", "1"]
    q7 <- counts["1", "1", "0"]
    q8 <- counts["1", "1", "1"]
    
    qhat <- c(q1, q2, q3, q4, q5, q6, q7, q8) / sum(counts)
    
    if ((qhat[3] + qhat[4]) > 0 && (qhat[7] + qhat[8]) > 0) {
      EY_r1_a0 <- qhat[4] / (qhat[3] + qhat[4])
      EY_r1_a1 <- qhat[8] / (qhat[7] + qhat[8])
      P_r1     <- qhat[3] + qhat[4] + qhat[7] + qhat[8]
      
      nb.ms[i] <- (EY_r1_a0 - EY_r1_a1 - t.grd[i]) * P_r1
    } else {
      nb.ms[i] <- NA
    }
    
    # Dirichlet posterior
    alpha <- c(q1, q2, q3, q4, q5, q6, q7, q8) + 1
    
    gamma_draws <- matrix(
      rgamma(m * 8, shape = rep(alpha, each = m)),
      nrow = m, ncol = 8, byrow = FALSE
    )
    
    dirichlet_draws <- gamma_draws / rowSums(gamma_draws)
    
    for (j in 1:m) {
      q <- dirichlet_draws[j, ]
      
      if ((q[3] + q[4]) > 0 && (q[7] + q[8]) > 0) {
        EY_r1_a0 <- q[4] / (q[3] + q[4])
        EY_r1_a1 <- q[8] / (q[7] + q[8])
        P_r1     <- q[3] + q[4] + q[7] + q[8]
        
        nb.posterior[j, i] <- (EY_r1_a0 - EY_r1_a1 - t.grd[i]) * P_r1
      } else {
        nb.posterior[j, i] <- NA
      }
    }
    
    # Posterior summaries
    nb.mean[i] <- mean(nb.posterior[, i], na.rm = TRUE)
    nb.low[i]  <- quantile(nb.posterior[, i], 0.025, na.rm = TRUE)
    nb.high[i] <- quantile(nb.posterior[, i], 0.975, na.rm = TRUE)
  }
  
  list(
    t.grd = t.grd,
    nb.plugin = nb.ms,
    nb.posterior = nb.posterior,
    nb.posterior.mean = nb.mean,
    nb.lower = nb.low,
    nb.upper = nb.high
  )
}

##### Calculating the bootstrap standard errors for the ATT Based Estimator ####
bootstrap_NB_ATTB <- function(y, a, pb, t.grd = NULL, B = 1000) {
  
  n <- length(y)
  
  if (is.null(t.grd)) {
    t.grd <- sort(unique(pb))
  }
  
  nb.boot <- matrix(NA, nrow = B, ncol = length(t.grd))
  
  for (b in 1:B) {
    
    idx <- sample(1:n, size = n, replace = TRUE)
    
    boot.out <- NB.ATTB(y = y[idx],
                      a = a[idx],
                      pb = pb[idx],
                      t.grd = t.grd)
    
    nb.boot[b, ] <- boot.out$nb.model
  }
  
  nb.lower <- apply(nb.boot, 2, quantile, probs = 0.025, na.rm = TRUE)
  nb.upper <- apply(nb.boot, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  list(
    t.grd = t.grd,
    nb.boot = nb.boot,
    
    ci.lower = nb.lower,
    ci.upper = nb.upper
  )
}

#Model-based NB. Only valid if pb=h. Does not work with y and a
NB.mb <- function(pb, t.grd=NULL) {
  
  ### doing quickly, but ...
  ### should check all three inputs are same length 
  ### should check pb is numeric, in [-1,1]
  
  ##  very fine grid
  if(is.null(t.grd))
  {
    t.grd <- sort(unique(pb))
  }
  
  nb.model <- rep(NA, length(t.grd))  ### will be output
  
  for (i in 1:length(t.grd)) {
    
    ### who is above threshold?
    ndx <- which(pb >= t.grd[i])
    nb.model[i] <- (mean(pb[ndx])-t.grd[i])*mean(pb > t.grd[i]) 
  }
  
  list(t.grd=t.grd, nb.model=nb.model, 
       nb.all=mean(pb)-t.grd)
}  

##### Lets now plot some curves ####

t.grd <- (1:100)/1000

cb <- NB.CB(TBP_eval$y, TBP_eval$a, TBP_eval$pb, t.grd)
attb <- NB.ATTB(TBP_eval$y, TBP_eval$a, TBP_eval$pb, t.grd)
mb <- NB.mb(TBP_eval$pb, t.grd)


bayes.cb.unshared <- bayesian_NB_CB_Unshared(TBP_eval$y, TBP_eval$a, TBP_eval$pb, t.grd)
bayes.attb.unshared <- bayesian_NB_ATTB_Unshared(TBP_eval$y, TBP_eval$a, TBP_eval$pb, t.grd)

se.cb <- compute_se_cb(TBP_eval$y, TBP_eval$a, TBP_eval$pb, t.grd)
se.attb <- compute_se_attb(TBP_eval$y, TBP_eval$a, TBP_eval$pb, t.grd)

bootstrapped_cb <- bootstrap_NB_CB(TBP_eval$y, TBP_eval$a, TBP_eval$pb, t.grd)
bootstrapped_attb <- bootstrap_NB_ATTB(TBP_eval$y, TBP_eval$a, TBP_eval$pb, t.grd)


### Time to decorate all our plots! 


plot(cb$t.grd, cb$nb.model, type='l', col='blue')
lines(attb$t.grd, attb$nb.model, col='green')
lines(mb$t.grd, mb$nb.model,col='black')

### what if we treat all 
lines(attb$t.grd, attb$nb.all, col='hotpink1')
lines(cb$t.grd, cb$nb.all,col= 'hotpink3')

lines(mb$t.grd, mb$nb.all,col= 'magenta4')


lines(se.cb$t.grd, se.cb$ci.lower, col = 'darkblue', lty = 2)
lines(se.cb$t.grd, se.cb$ci.upper, col = 'darkblue', lty = 2)


lines(se.attb$t.grd, se.attb$ci.lower, col = 'darkgreen', lty = 2)
lines(se.attb$t.grd, se.attb$ci.upper, col = 'darkgreen', lty = 2)


lines(bayes.cb.unshared$t.grd, bayes.cb.unshared$nb.posterior.mean, col='darkorange')
lines(bayes.cb.unshared$t.grd, bayes.cb.unshared$ci.lower, col='darkorange3', lty = 2)
lines(bayes.cb.unshared$t.grd, bayes.cb.unshared$ci.upper, col='darkorange3',lty = 2)


lines(bayes.attb.unshared$t.grd, bayes.attb.unshared$nb.posterior.mean, col='green3')
lines(bayes.attb.unshared$t.grd, bayes.attb.unshared$nb.lower, col='darkgreen3', lty = 2)
lines(bayes.attb.unshared$t.grd, bayes.attb.unshared$nb.upper, col='darkgreen3', lty = 2)


lines(bootstrapped_cb$t.grd, bootstrapped_cb$ci.lower, col = 'magenta')
lines(bootstrapped_cb$t.grd, bootstrapped_cb$ci.upper, col = 'magenta')

lines(bootstrapped_attb$t.grd, bootstrapped_attb$ci.lower, col = 'darkorchid')
lines(bootstrapped_attb$t.grd, bootstrapped_attb$ci.upper, col = 'darkorchid')


legend("topright",
       legend = c(
         "PG", 
         "PG 95% CI (asymptotic)", 
         "PG posterior mean", 
         "PG 95% CI (Bayes)", 
         "PG 95% Bootstrap Interval",
         "MS", 
         "MS 95% CI (asymptotic)", 
         "MS posterior mean", 
         "MS 95% CI (Bayes)", 
         "MS 95% Bootstrap Interval",
         "MB"
       ),
       col = c(
         "blue",          # CB plugin
         "darkblue",      # CB asymptotic CI
         "darkorange",    # CB posterior mean
         "darkorange3",   # PCB posterior CI
         "magenta",       # CB Bootstrap Interval
         "green",         # ATTB plugin
         "darkgreen",     # ATTB asymptotic CI
         "green",         # ATTB posterior mean
         "darkgreen",     # ATTB posterior CI
         "seagreen",      # ATTB Bootstrap Interval 
         "black"          # ATTB plugin
       ),
       lty = c(
         1, 2, 1, 2, 1, 2, 1, 2, 1
       ),
       lwd = c(
         1, 1, 2, 1, 1, 1, 2, 1, 1
       ),
       title = "Estimator Type",
       cex = 0.8,
       bty = "n"
)

##### Chatgpt is doing something?? 

library(dplyr)
library(tidyr)
library(ggplot2)

# Plugin estimates
df_plugin <- bind_rows(
  data.frame(t = cb$t.grd, nb = cb$nb.model, estimator = "CB", method = "Plugin"),
  data.frame(t = attb$t.grd, nb = attb$nb.model, estimator = "ATTB", method = "Plugin"),
  data.frame(t = mb$t.grd, nb = mb$nb.model, estimator = "MB", method = "Plugin")
)

# Treat All lines (flat)
df_treatall <- bind_rows(
  data.frame(t = cb$t.grd, nb = cb$nb.all, estimator = "CB", method = "Treat All"),
  data.frame(t = attb$t.grd, nb = attb$nb.all, estimator = "ATTB", method = "Treat All"),
  data.frame(t = mb$t.grd, nb = mb$nb.all, estimator = "MB", method = "Treat All")
)

# Asymptotic CIs
df_asymp <- bind_rows(
  data.frame(t = se.cb$t.grd, nb = se.cb$ci.lower, estimator = "CB", method = "Asymptotic Lower"),
  data.frame(t = se.cb$t.grd, nb = se.cb$ci.upper, estimator = "CB", method = "Asymptotic Upper"),
  data.frame(t = se.attb$t.grd, nb = se.attb$ci.lower, estimator = "ATTB", method = "Asymptotic Lower"),
  data.frame(t = se.attb$t.grd, nb = se.attb$ci.upper, estimator = "ATTB", method = "Asymptotic Upper")
)

# Bayesian posterior means and intervals
df_bayes <- bind_rows(
  data.frame(t = bayes.cb.unshared$t.grd, nb = bayes.cb.unshared$nb.posterior.mean, estimator = "CB", method = "Bayes Mean"),
  data.frame(t = bayes.cb.unshared$t.grd, nb = bayes.cb.unshared$ci.lower, estimator = "CB", method = "Bayes Lower"),
  data.frame(t = bayes.cb.unshared$t.grd, nb = bayes.cb.unshared$ci.upper, estimator = "CB", method = "Bayes Upper"),
  
  data.frame(t = bayes.attb.unshared$t.grd, nb = bayes.attb.unshared$nb.posterior.mean, estimator = "ATTB", method = "Bayes Mean"),
  data.frame(t = bayes.attb.unshared$t.grd, nb = bayes.attb.unshared$nb.lower, estimator = "ATTB", method = "Bayes Lower"),
  data.frame(t = bayes.attb.unshared$t.grd, nb = bayes.attb.unshared$nb.upper, estimator = "ATTB", method = "Bayes Upper")
)

# Bootstrap intervals
df_boot <- bind_rows(
  data.frame(t = bootstrapped_cb$t.grd, nb = bootstrapped_cb$ci.lower, estimator = "CB", method = "Bootstrap Lower"),
  data.frame(t = bootstrapped_cb$t.grd, nb = bootstrapped_cb$ci.upper, estimator = "CB", method = "Bootstrap Upper"),
  
  data.frame(t = bootstrapped_attb$t.grd, nb = bootstrapped_attb$ci.lower, estimator = "ATTB", method = "Bootstrap Lower"),
  data.frame(t = bootstrapped_attb$t.grd, nb = bootstrapped_attb$ci.upper, estimator = "ATTB", method = "Bootstrap Upper")
)

# Combine all
df_all <- bind_rows(df_plugin, df_treatall, df_asymp, df_bayes, df_boot)


ggplot(df_all, aes(x = t, y = nb, color = method, linetype = method)) +
  geom_line(size = 0.5) +
  facet_wrap(~ estimator, nrow = 1) +
  labs(
    title = "Net Benefit Curves by Estimator and Method",
    x = "Threshold (T)",
    y = "Net Benefit",
    color = "Method",
    linetype = "Method"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c(
    "Plugin" = "black",
    "Treat All" = "blue",
    "Asymptotic Lower" = "dodgerblue4",
    "Asymptotic Upper" = "dodgerblue4",
    "Bayes Mean" = "hotpink",
    "Bayes Lower" = "hotpink4",
    "Bayes Upper" = "hotpink4",
    "Bootstrap Lower" = "darkorchid4",
    "Bootstrap Upper" = "darkorchid4"
  )) +
  scale_linetype_manual(values = c(
    "Plugin" = "solid",
    "Treat All" = "solid",
    "Asymptotic Lower" = "dashed",
    "Asymptotic Upper" = "dashed",
    "Bayes Mean" = "solid",
    "Bayes Lower" = "dashed",
    "Bayes Upper" = "dashed",
    "Bootstrap Lower" = "dashed",
    "Bootstrap Upper" = "dashed"
  )) +
  coord_cartesian(ylim = c(-0.005, NA))











