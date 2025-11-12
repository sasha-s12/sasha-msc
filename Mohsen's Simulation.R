
### write about survival outcome and confouding in discussion maybe

### Try it on real data -> use simulation study for general cases

## needs to be aware of randomization ratio 

## try the new bayesian thing 

## mse coverage probability <- try for both confidence interval formulas looking for same average length and coverage probabiltiy

set.seed(123)

#### Treating Probability A is known as a fixed parameter
gen_data <- function(n=1000, r=0.25)
{
  x <- rnorm(n,1,1) 
  a <- rbinom(n,1,r) 
  p0 <- 1/(1+exp(-(-2+x)))
  p1 <- 1/(1+exp(-(-2+x-0.5)))
  Y0 <- rbinom(n, 1, p0)
  Y1 <- rbinom(n, 1, p1)
  Y <- a*Y1+(1-a)*Y0
  h <- p0-p1
  df <- data.frame(h=h, a=a, Y=Y, p0=p0)
  df <- df[order(df$h),]
  df
}

#### Paul's Estimator ####
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

bayesian_NB_PG_Shared <- function(y, a, pb, t.grd, m = 1000) {
  
  n <- length(y)
  nb.posterior <- matrix(NA, nrow = m, ncol = length(t.grd))
  nb.pg <- numeric(length(t.grd))
  nb.mean <- numeric(length(t.grd))
  nb.lower <- numeric(length(t.grd))
  nb.upper <- numeric(length(t.grd))
  
  # Shared Uniform(0,1) draws to be used for inverse CDF gamma transform
  unif_base <- matrix(runif(m * 8), nrow = m, ncol = 8)
  
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
    
    # Transform shared uniform draws to gamma using threshold-specific alpha
    gamma_draws <- matrix(NA, nrow = m, ncol = 8)
    for (k in 1:8) {
      gamma_draws[, k] <- qgamma(unif_base[, k], shape = alpha[k], scale = 1)
    }
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




#### ATT Based Estimator ####
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

#### Calculating the Bayesian Standard Errors for the ATT-Based Estimator #### 
bayesian_NB_MS_Shared <- function(y, a, pb, t.grd, m = 1000) {
  
  nb.posterior <- matrix(NA, nrow = m, ncol = length(t.grd))
  nb.plugin <- numeric(length(t.grd))      # Plug-in MS estimator
  nb.mean   <- numeric(length(t.grd))      # Posterior mean
  nb.low    <- numeric(length(t.grd))      # 2.5% quantile
  nb.high   <- numeric(length(t.grd))      # 97.5% quantile
  
  # Step 1: Shared uniform(0,1) draws
  unif_base <- matrix(runif(m * 8), nrow = m, ncol = 8)
  
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
    
    # --- Plug-in MS estimator ---
    qhat <- as.numeric(c(q1, q2, q3, q4, q5, q6, q7, q8)) / sum(counts)
    
    denom_a0 <- qhat[3] + qhat[4]
    denom_a1 <- qhat[7] + qhat[8]
    
    if (denom_a0 > 0 && denom_a1 > 0) {
      EY_r1_a0 <- qhat[4] / denom_a0
      EY_r1_a1 <- qhat[8] / denom_a1
      P_r1     <- denom_a0 + denom_a1
      
      nb.plugin[i] <- (EY_r1_a0 - EY_r1_a1 - t.grd[i]) * P_r1
    } else {
      nb.plugin[i] <- NA
    }
    
    # --- Dirichlet posterior ---
    alpha <- c(q1, q2, q3, q4, q5, q6, q7, q8) + 1
    
    # Transform shared uniforms to gamma
    gamma_draws <- matrix(NA, nrow = m, ncol = 8)
    for (k in 1:8) {
      gamma_draws[, k] <- qgamma(unif_base[, k], shape = alpha[k], scale = 1)
    }
    
    dirichlet_draws <- gamma_draws / rowSums(gamma_draws)
    
    for (j in 1:m) {
      q <- dirichlet_draws[j, ]
      
      denom_a0 <- q[3] + q[4]
      denom_a1 <- q[7] + q[8]
      
      if (denom_a0 > 0 && denom_a1 > 0) {
        EY_r1_a0 <- q[4] / denom_a0
        EY_r1_a1 <- q[8] / denom_a1
        P_r1     <- denom_a0 + denom_a1
        
        nb.posterior[j, i] <- (EY_r1_a0 - EY_r1_a1 - t.grd[i]) * P_r1
      }
    }
    
    nb.mean[i] <- mean(nb.posterior[, i], na.rm = TRUE)
    nb.low[i]  <- quantile(nb.posterior[, i], 0.025, na.rm = TRUE)
    nb.high[i] <- quantile(nb.posterior[, i], 0.975, na.rm = TRUE)
  }
  
  list(
    t.grd = t.grd,
    nb.plugin = nb.plugin,
    nb.posterior = nb.posterior,
    nb.posterior.mean = nb.mean,
    nb.lower = nb.low,
    nb.upper = nb.high
  )
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


df <- gen_data(10^5)

t.grd <- (1:100)/1000
pg <- NB.PG(df$Y, df$a, df$h, t.grd)
ms <- NB.MS(df$Y, df$a, df$h, t.grd)
mb <- NB.mb(df$h, t.grd)

plot(pg$t.grd, pg$nb.model, type='l', col='blue')
lines(ms$t.grd, ms$nb.model, type = 'l', col='green')
lines(mb$t.grd, mb$nb.model,col='black')

#legend("topright",legend = c("PG", "MS", "MB"),col = c("blue", "green", "black"), lty = 1,  title = "Estimator")


bayes.pg.shared <- bayesian_NB_PG_Shared(df$Y, df$a, df$h, t.grd)

bayes.ms.shared <- bayesian_NB_MS_Shared(df$Y, df$a, df$h, t.grd)


se.pg <- compute_se_pg(df$Y, df$a, df$h, t.grd)
se.ms <- compute_se_ms(df$Y, df$a, df$h, t.grd)


lines(se.pg$t.grd, se.pg$ci.lower, col = 'darkblue', lty = 2)
lines(se.pg$t.grd, se.pg$ci.upper, col = 'darkblue', lty = 2)


lines(se.ms$t.grd, se.ms$ci.lower, col = 'darkgreen', lty = 2)
lines(se.ms$t.grd, se.ms$ci.upper, col = 'darkgreen', lty = 2)


#plot(bayes.ms.unshared$t.grd, bayes.ms.unshared$nb.plugin, type='l', col='blue')

#lines(bayes.ms.unshared$t.grd, bayes.ms.unshared$nb.plugin, type='l', col="darkred")

lines(bayes.ms.unshared$t.grd, bayes.ms.unshared$nb.posterior.mean, col='green')

lines(bayes.ms.unshared$t.grd, bayes.ms.unshared$nb.lower, col='darkgreen', lty = 2)

lines(bayes.ms.unshared$t.grd, bayes.ms.unshared$nb.upper, col='darkgreen', lty = 2)


#lines(bayes.pg.unshared$t.grd, bayes.pg.unshared$nb.plugin, type='l', col="blue")

lines(bayes.pg.unshared$t.grd, bayes.pg.unshared$nb.posterior.mean, col='darkorange')

lines(bayes.pg.unshared$t.grd, bayes.pg.unshared$ci.lower, col='darkorange3', lty = 2)

lines(bayes.pg.unshared$t.grd, bayes.pg.unshared$ci.upper, col='darkorange3',lty = 2)


bootstrapped_ms <- bootstrap_NB_MS(df$Y, df$a, df$h, t.grd)

bootstrapped_pg <- bootstrap_NB_PG(df$Y, df$a, df$h, t.grd)

lines(bootstrapped_pg$t.grd, bootstrapped_pg$ci.lower, col = 'magenta')
lines(bootstrapped_pg$t.grd, bootstrapped_pg$ci.upper, col = 'magenta')

lines(bootstrapped_ms$t.grd, bootstrapped_ms$ci.lower, col = 'seagreen')
lines(bootstrapped_ms$t.grd, bootstrapped_ms$ci.upper, col = 'seagreen')


lines(bayes.pg.shared$t.grd, bayes.pg.shared$nb.posterior.mean, col='red')

lines(bayes.pg.shared$t.grd, bayes.pg.shared$ci.lower, col='plum3', lty = 2)

lines(bayes.pg.shared$t.grd, bayes.pg.shared$ci.upper, col='plum3',lty = 2)

lines(bayes.pg.shared$t.grd, bayes.pg.shared$nb.plugin, col='red',lty = 2)



lines(bayes.ms.shared$t.grd, bayes.ms.shared$nb.posterior.mean, col='red')

lines(bayes.ms.shared$t.grd, bayes.ms.shared$ci.lower, col='plum3', lty = 2)

lines(bayes.ms.shared$t.grd, bayes.ms.shared$ci.upper, col='plum3',lty = 2)

lines(bayes.ms.shared$t.grd, bayes.ms.shared$nb.plugin, col='red',lty = 2)


plot(bayes.pg.shared$t.grd, bayes.pg.shared$nb.posterior.mean, type='l', col='blue')

lines(pg$t.grd, ms$nb.model, col='pink')
lines(ms$t.grd, ms$nb.model, col='green')
lines(mb$t.grd, mb$nb.model,col='black')


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
         "blue",          # PG plugin
         "darkblue",      # PG asymptotic CI
         "darkorange",    # PG posterior mean
         "darkorange3",   # PG posterior CI
         "magenta",       # PG Bootstrap Interval
         "green",         # MS plugin
         "darkgreen",     # MS asymptotic CI
         "green",         # MS posterior mean
         "darkgreen",     # MS posterior CI
         "seagreen",      # MS Bootstrap Interval 
         "black"          # MB plugin
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




