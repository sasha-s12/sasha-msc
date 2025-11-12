
set.seed(123)

library(dplyr)

##### Estimators #####

#### Paul's Estimator
NB.PG <- function(y, a, pb, t.grd=NULL) {
  
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

#### Mohsen's Estimator
NB.MS <- function(y, a, pb, t.grd=NULL) {
  
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

#Model-based estimator (makes sense only for when model is correct)
NB.param <- function(y, a, pb, t.grd=NULL) {
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
       nb.all=mean(pb-t.grd))
}





t.grd <- (1:100) / 1000





##### CI Calculators ####

### Calculating the Standard Errors for the Congruent Based Estimator ####
compute_se_pg <- function(y, a, pb, t.grd) {
  
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
bayesian_NB_PG_Unshared <- function(y, a, pb, t.grd, m = 1000) {
  
  n <- length(y)
  nb.posterior <- matrix(NA, nrow = m, ncol = length(t.grd))
  nb.pg <- numeric(length(t.grd))     # plug-in
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
    
    nb.pg[i] <- benefit - (harm_ctrl + harm_trt) - t.grd[i] * p
    
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
    nb.plugin = nb.pg,
    nb.posterior = nb.posterior,
    nb.posterior.mean = nb.mean,
    ci.lower = nb.lower,
    ci.upper = nb.upper
  )
}


##### Calculating the bootstrap standard errors for the congruent based estimator ####
bootstrap_NB_PG <- function(y, a, pb, t.grd = NULL, B = 1000) {
  
  n <- length(y)
  
  if (is.null(t.grd)) {
    t.grd <- sort(unique(pb))
  }
  
  nb.boot <- matrix(NA, nrow = B, ncol = length(t.grd))
  
  for (b in 1:B) {
    
    idx <- sample(1:n, size = n, replace = TRUE)
    
    boot.out <- NB.PG(y = y[idx],
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

compute_se_ms <- function(y, a, pb, t.grd) {
  
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

#### Calculating the Bayesian Standard Errors for the ATT-Based Estimator Unshared ####
bayesian_NB_MS_Unshared <- function(y, a, pb, t.grd, m = 1000) {
  
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
bootstrap_NB_MS <- function(y, a, pb, t.grd = NULL, B = 1000) {
  
  n <- length(y)
  
  if (is.null(t.grd)) {
    t.grd <- sort(unique(pb))
  }
  
  nb.boot <- matrix(NA, nrow = B, ncol = length(t.grd))
  
  for (b in 1:B) {
    
    idx <- sample(1:n, size = n, replace = TRUE)
    
    boot.out <- NB.MS(y = y[idx],
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



##### Code to calculate the true value of H using numerical integration ####

# Code to generate the "true" data set using a randomization ratio of 1:3 
gen_data_true <- function(n=1000000, r=1/4)
{
  x <- rnorm(n,1,1) 
  a <- rbinom(n,1,r) 
  p0 <- 1/(1+exp(-(-2+x)))
  p1 <- 1/(1+exp(-(-2+x-0.5)))
  Y0 <- rbinom(n, 1, p0)
  Y1 <- rbinom(n, 1, p1)
  Y <- a*Y1+(1-a)*Y0
  h <- (p0-p1)
  df <- data.frame(h=h, a=a, Y=Y, p0=p0, x=x)
  df <- df[order(df$h),]
  df
}

df <- gen_data_true()

simulate_one_true <- function(t.grd) {
  
  df <- gen_data_true()
  
  pg <- NB.PG(df$Y, df$a, df$h, t.grd)
  ms <- NB.MS(df$Y, df$a, df$h, t.grd)
  mb <- NB.mb(df$h, t.grd)
  
  data.frame(
    threshold = t.grd,
    NB_PG = pg$nb.model,
    NB_MS = ms$nb.model,
    NB_MB = mb$nb.model
  )
}

NB_true <- simulate_one_true(t.grd)

NB_true


##### Now we are creating the samples ####

### Code to generate a sample dataset ###
gen_data <- function(n=100000, r=1/4)
{
  x <- rnorm(n,1,1) 
  a <- rbinom(n,1,r) 
  p0 <- 1/(1+exp(-(-2+x)))
  p1 <- 1/(1+exp(-(-2+x-0.5)))
  Y0 <- rbinom(n, 1, p0)
  Y1 <- rbinom(n, 1, p1)
  Y <- a*Y1+(1-a)*Y0
  h <- (p0-p1)
  df <- data.frame(h=h, a=a, Y=Y, p0=p0, x=x)
  df <- df[order(df$h),]
  df
}

simulate_one <- function(t.grd) {
  
  df <- gen_data(1000)

  pg <- NB.PG(df$Y, df$a, df$h, t.grd)
  ms <- NB.MS(df$Y, df$a, df$h, t.grd)
  mb <- NB.mb(df$h, t.grd)
  
  data.frame(
    threshold = t.grd,
    NB_PG = pg$nb.model,
    NB_MS = ms$nb.model,
    NB_MB = mb$nb.model
  )
}

run_simulation <- function(n_sim) {
  
  results <- lapply(1:n_sim, function(i) {
    
    sim_result <- simulate_one(t.grd)
    sim_result$sim <- i
    sim_result }) %>% dplyr::bind_rows()
  
  return(results)
}

sim_results <- run_simulation(n_sim = 500)

sim_results[sim_results$sim == 500, ]

##### Now to calculate the MSEs ####

# First we will merge the true dataset and the simulated dataset 

merged <- merge(sim_results, NB_true, by = "threshold", suffixes = c("", "_true"))

# Now we are calculating the MSE per threshold per estimator 

merged$SE_PG <- (merged$NB_PG - merged$NB_PG_true)^2
merged$SE_MS <- (merged$NB_MS - merged$NB_MS_true)^2
merged$SE_MB <- (merged$NB_MB - merged$NB_MB_true)^2


### Now we take the average SE over each threshold

mse_PG <- tapply(merged$SE_PG, merged$threshold, mean)
mse_MS <- tapply(merged$SE_MS, merged$threshold, mean)
mse_MB <- tapply(merged$SE_MB, merged$threshold, mean)

mse_summary <- data.frame(
  threshold = as.numeric(names(mse_PG)),
  MSE_PG = mse_PG,
  MSE_MS = mse_MS,
  MSE_MB = mse_MB
)


##### Making some plots for the MSE at each threshold comparing the three curves####

plot(mse_summary$threshold, mse_summary$MSE_PG,
     type = "l", col = "red", lwd = 2,
     xlab = "Threshold", ylab = "MSE",
     main = "MSE of Net Benefit Estimators")  

lines(mse_summary$threshold, mse_summary$MSE_MS, col = "blue", lwd = 2)
lines(mse_summary$threshold, mse_summary$MSE_MB, col = "green", lwd = 2)

legend("topright",                     # Legend position
       legend = c("PG", "MS", "MB"),   # Labels for each line
       col = c("red", "blue", "green"),# Colors used in the plot
       lty = 1,                        # Line type (1 = solid)
       lwd = 2)                        # Line width


plot(mse_summary$threshold, mse_summary$MSE_MB, type = "l", col = "green", lwd = 2)


#### Comparing the NB for the noisy one vs the true one across thresholds only do this after we change the h value manualy  

# random data set
NB_noisy <- sim_results[sim_results$sim == 1, ]

# Merge with true NB
nb_compare <- merge(NB_noisy, NB_true, by = "threshold", suffixes = c("_noisy", "_true"))

# Plot all 6 curves correctly
plot(nb_compare$threshold, nb_compare$NB_PG_true, type = "l", col = "red", lwd = 2,
     ylim = range(nb_compare$NB_PG_true,
                  nb_compare$NB_PG_noisy,
                  nb_compare$NB_MS_true,
                  nb_compare$NB_MS_noisy,
                  nb_compare$NB_MB_true,
                  nb_compare$NB_MB_noisy,
                  na.rm = TRUE),
     xlab = "Threshold", ylab = "Net Benefit",
     main = "True vs Noisy NB Curves (PG, MS, MB)")


lines(nb_compare$threshold, nb_compare$NB_PG_noisy, col = "red", lty = 2, lwd = 2)
lines(nb_compare$threshold, nb_compare$NB_MS_true, col = "blue", lwd = 2)
lines(nb_compare$threshold, nb_compare$NB_MS_noisy, col = "blue", lty = 2, lwd = 2)
lines(nb_compare$threshold, nb_compare$NB_MB_true, col = "green", lwd = 2)
lines(nb_compare$threshold, nb_compare$NB_MB_noisy, col = "green", lty = 2, lwd = 2)

legend("topright",
       legend = c("PG True", "PG Noisy", "MS True", "MS Noisy", "MB True", "MB Noisy"),
       col = c("red", "red", "blue", "blue", "green", "green"),
       lty = c(1, 2, 1, 2, 1, 2),
       lwd = 2)






##### Calculating coverage probabilities and confidence interval length ####

evaluate_coverage <- function(n_sim, t.grd, true_vals, ci_function, estimator_label, method_label) {
  
  coverage_results <- matrix(0, nrow = length(t.grd), ncol = n_sim)
  ci_lengths <- matrix(NA, nrow = length(t.grd), ncol = n_sim)
  
  for (i in 1:n_sim) {
    
    df <- gen_data(1000)
    y <- df$Y
    a <- df$a
    pb <- df$h  # Assuming using h as the predicted benefit
    
    ci <- ci_function(y, a, pb, t.grd)
    
    # True NB vector
    theta_true <- true_vals$NB_true
    
    # Whether the CI contains the true value
    coverage_results[, i] <- (theta_true >= ci$ci.lower) & (theta_true <= ci$ci.upper)
    
    # CI length
    ci_lengths[, i] <- ci$ci.upper - ci$ci.lower
  }
  
  coverage <- rowMeans(coverage_results, na.rm = TRUE)
  avg_length <- rowMeans(ci_lengths, na.rm = TRUE)
  
  data.frame(
    threshold = t.grd,
    coverage = coverage,
    avg_ci_length = avg_length,
    method = method_label,
    estimator = estimator_label
  )
}

evaluate_coverage_parallel <- function(n_sim, t.grd, true_vals, ci_function, 
                                       estimator_label, method_label, n_cores = detectCores() - 1) {
  
  # wrapper function for one simulation
  single_sim <- function(i) {
    df <- gen_data(1000)
    y <- df$Y
    a <- df$a
    pb <- df$h
    
    ci <- ci_function(y, a, pb, t.grd)
    
    coverage <- (true_vals$NB_true >= ci$ci.lower) & (true_vals$NB_true <= ci$ci.upper)
    ci_len <- ci$ci.upper - ci$ci.lower
    
    list(sim_id = i, coverage = coverage, ci_length = ci_len)
  }
  
  results <- mclapply(1:n_sim, single_sim, mc.cores = n_cores) # returns nsim lists with scalar sim_id = i, two vectors of length t.grd
  
  coverage_matrix <- do.call(cbind, lapply(results, `[[`, "coverage")) ## the lapply function just extracts the coverage from each nsim list then cbinds it based on nsim
  ci_length_matrix <- do.call(cbind, lapply(results, `[[`, "ci_length")) 
  
  data.frame(
    threshold = t.grd,
    coverage = rowMeans(coverage_matrix, na.rm = TRUE),
    avg_ci_length = rowMeans(ci_length_matrix, na.rm = TRUE),
    method = method_label,
    estimator = estimator_label
  )
}


true_pg_df <- data.frame(threshold = NB_true$threshold, NB_true = NB_true$NB_PG)

pg_asymptotic <- evaluate_coverage_parallel(
  n_sim = 500,
  t.grd = t.grd,
  true_vals = true_pg_df,
  ci_function = compute_se_pg,       # this is the same function you used before
  estimator_label = "PG",
  method_label = "Asymptotic"
)

pg_bayes <- evaluate_coverage_parallel(
  500,                      # Number of simulations
  t.grd,                    # Threshold grid
  true_pg_df,              # True NB values for PG
  function(y, a, pb, t.grd) {
    out <- bayesian_NB_PG_Unshared(y, a, pb, t.grd)
    list(
      theta_hat = rep(NA, length(t.grd)),   # Not needed for coverage here
      ci.lower = out$ci.lower,
      ci.upper = out$ci.upper
    )
  },
  "PG",                    # Estimator label
  "Bayesian"              # CI method label
)

pg_bootstrap <- evaluate_coverage_parallel(
  n_sim = 500,
  t.grd = t.grd,
  true_vals = true_pg_df,
  ci_function = function(y, a, pb, t.grd) {
    out <- bootstrap_NB_PG(y, a, pb, t.grd, B = 100)  # Temporarily reduce B while testing
    list(theta_hat = rep(NA, length(t.grd)), ci.lower = out$ci.lower, ci.upper = out$ci.upper)
  },
  estimator_label = "PG",
  method_label = "Bootstrap"
)


pg_combined <- bind_rows(pg_asymptotic, pg_bootstrap, pg_bayes)


library(ggplot2)

ggplot(pg_combined, aes(x = threshold, y = coverage, color = method)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  ylim(0.9, 1) +
  labs(
    title = "Coverage Probability (C-B Estimator)",
    x = "Threshold",
    y = "Coverage Probability",
    color = "Method"    
  ) +
  theme_minimal()


ggplot(pg_combined, aes(x = threshold, y = avg_ci_length, color = method)) +
  geom_line(size = 1.2) +
  labs(title = "Average CI Length - C-B Estimator", x = "Threshold", y = "Avg CI Length",
       color = "Method"
       ) +
  theme_minimal()


true_ms_df <- data.frame(threshold = NB_true$threshold, NB_true = NB_true$NB_MS)

# Asymptotic for MS
ms_asymptotic <- evaluate_coverage_parallel(
  n_sim = 500,
  t.grd = t.grd,
  true_vals = true_ms_df,
  ci_function = compute_se_ms,
  estimator_label = "MS",
  method_label = "Asymptotic"
)

# Bayesian for MS
ms_bayes <- evaluate_coverage_parallel(
  500,
  t.grd,
  true_ms_df,
  function(y, a, pb, t.grd) {
    out <- bayesian_NB_MS_Unshared(y, a, pb, t.grd)
    list(
      theta_hat = rep(NA, length(t.grd)),
      ci.lower = out$nb.lower,
      ci.upper = out$nb.upper
    )
  },
  estimator_label = "MS",
  method_label = "Bayesian"
)

# Bootstrap for MS
ms_bootstrap <- evaluate_coverage_parallel(
  n_sim = 500,
  t.grd = t.grd,
  true_vals = true_ms_df,
  ci_function = function(y, a, pb, t.grd) {
    out <- bootstrap_NB_MS(y, a, pb, t.grd, B = 100)  # Adjust B as needed
    list(
      theta_hat = rep(NA, length(t.grd)),
      ci.lower = out$ci.lower,
      ci.upper = out$ci.upper
    )
  },
  estimator_label = "MS",
  method_label = "Bootstrap"
)

ms_combined <- bind_rows(ms_asymptotic, ms_bootstrap, ms_bayes)

ggplot(ms_combined, aes(x = threshold, y = coverage, color = method)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  ylim(0.9, 1) +
  labs(
    title = "Coverage Probability (ATT-B Estimator)",
    x = "Threshold",
    y = "Coverage Probability",
    color = "Method"    
  ) +
  theme_minimal()


ggplot(ms_combined, aes(x = threshold, y = avg_ci_length, color = method)) +
  geom_line(size = 1.2) +
  labs(title = "Average CI Length - ATT-B Estimator", x = "Threshold", y = "Avg CI Length",
       color = "Method"
  ) +
  theme_minimal()







