
library(dplyr)
library(parallel)

set.seed(123)

##### Treatment Thresholds 
t.grd <- (1:100) / 1000

##### Estimators #####

#### Congruent Based Estimator
NB.CG <- function(y, a, pb, t.grd=NULL) {
  
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

#### ATT Based Estimator
NB.ATT <- function(y, a, pb, t.grd=NULL) {
  
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

#### Model-based NB. Only valid if pb=h. Does not work with y and a
NB.MB <- function(pb, t.grd=NULL) {
  
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

#### Vickers-Style Estimator 
NB.VS <- function(y, a, pb, t.grd=NULL) {
  
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
    ### who is compliant at this threshold?
    
    ndx <- as.logical(a) == (pb > t.grd[i]) ### returns a list of the congruent individuals 
    
    nb.model[i] <- y.bar.cntrl - mean(y[ndx]) - t.grd[i]*mean(pb > t.grd[i])
    
  }
  
  list(t.grd=t.grd, nb.model=nb.model, 
       nb.all=y.bar.cntrl-y.bar.trt-t.grd)  
} 

##### CI Calculators 

##### Congruent Based Estimators ####

### Calculating the Standard Errors for the Congruent Based Estimator Delta Method 
delta_SE_CG <- function(y, a, pb, t.grd) {
  
  n <- length(y)
  se.cg <- numeric(length(t.grd))
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
    se.cg[i] <- sqrt(as.numeric(var_theta))
    
    ci.lower[i] <- theta_hat[i] - 1.96 * se.cg[i]
    ci.upper[i] <- theta_hat[i] + 1.96 * se.cg[i]
  }
  
  return(list(
    t.grd = t.grd,
    theta_hat = theta_hat,
    se = se.cg,
    ci.lower = ci.lower,
    ci.upper = ci.upper
  ))
}

### Calculating the bootstrap standard errors for the congruent based estimator 
bootstrap_SE_CG <- function(y, a, pb, t.grd = NULL, B = 1000) {
  
  n <- length(y)
  
  if (is.null(t.grd)) {
    t.grd <- sort(unique(pb))
  }
  
  nb.boot <- matrix(NA, nrow = B, ncol = length(t.grd))
  
  for (b in 1:B) {
    
    idx <- sample(1:n, size = n, replace = TRUE)
    
    boot.out <- NB.CG(y = y[idx],
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

### Calculating the Bayesian posterior distribution of net benefit
bayesian_CG_Shared <- function(y, a, pb, t.grd, m = 1000) {
  
  n <- length(y)
  nb.posterior <- matrix(NA, nrow = m, ncol = length(t.grd))
  nb.cg <- numeric(length(t.grd))
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
    
    nb.cg[i] <- benefit - (harm_ctrl + harm_trt) - t.grd[i] * p
    
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
    nb.plugin = nb.cg,
    nb.posterior = nb.posterior,
    nb.posterior.mean = nb.mean,
    ci.lower = nb.lower,
    ci.upper = nb.upper
  )
}


##### ATT Based Estimators ####
### Calculating the Delta Method Standard Errors for the ATT Based Estimator 
delta_SE_ATT <- function(y, a, pb, t.grd) {
  
  n <- length(y)
  se.att <- numeric(length(t.grd))
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
    se.att[i] <- sqrt(as.numeric(var_theta))
    
    # Confidence intervals
    ci.lower[i] <- theta_hat[i] - 1.96 * se.att[i]
    ci.upper[i] <- theta_hat[i] + 1.96 * se.att[i]
  }
  
  return(list(
    t.grd = t.grd,
    theta_hat = theta_hat,
    se = se.att,
    ci.lower = ci.lower,
    ci.upper = ci.upper
  ))
}

### Calculating the bootstrap standard errors for the ATT Based Estimator 
bootstrap_SE_ATT <- function(y, a, pb, t.grd = NULL, B = 1000) {
  
  n <- length(y)
  
  if (is.null(t.grd)) {
    t.grd <- sort(unique(pb))
  }
  
  nb.boot <- matrix(NA, nrow = B, ncol = length(t.grd))
  
  for (b in 1:B) {
    
    idx <- sample(1:n, size = n, replace = TRUE)
    
    boot.out <- NB.ATT(y = y[idx],
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

### Calculating the bootstrap standard errors for the ATT Based Estimator 
bayesian_ATT_Shared <- function(y, a, pb, t.grd, m = 1000) {
  
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


##### Code to calculate the true value of H using numerical integration ####

# Simulating the TRUE RCT using a randomization ratio of r # 
gen_data <- function(n, r, z = 1)
{
  x <- rnorm(n,1,1) 
  a <- rbinom(n,1,r) 
  p0 <- 1/(1+exp(-(-2+x)))
  p1 <- 1/(1+exp(-(-2+x-0.5)))
  Y0 <- rbinom(n, 1, p0)
  Y1 <- rbinom(n, 1, p1)
  Y <- a*Y1+(1-a)*Y0
  h <- (p0-p1)/z # M True value of H 
  df <- data.frame(h=h, a=a, Y=Y, p0=p0, x=x)
  df <- df[order(df$h),]
  df
}

# Simulating the true NB using a dataset # 
# Simulating the true NB using a dataset
simulate_one <- function(n, t.grd, r, z = 1, vs = FALSE) {
  
  df <- gen_data(n = n, r = r, z = z)
  
  cg  <- NB.CG(df$Y, df$a, df$h, t.grd)
  att <- NB.ATT(df$Y, df$a, df$h, t.grd)
  mb  <- NB.MB(df$h, t.grd)
  
  # base data frame
  out <- data.frame(
    threshold = t.grd,
    NB_CG = cg$nb.model,
    NB_ATT = att$nb.model,
    NB_MB = mb$nb.model
  )
  
  # conditionally add NB_VS
  if (vs == TRUE) {
    vs_out <- NB.VS(df$Y, df$a, df$h, t.grd)
    out$NB_VS <- vs_out$nb.model
  }
  
  out 
}

 NB_true_quarter <- simulate_one(n = 1000000, t.grd, r = 1/4, vs = TRUE) # 1:3 randomization 

##### Time to calculate my best friend Coverage probabilities and CI Length  ####

evaluate_coverage <- function(n_sim, t.grd, true_vals, ci_function, 
                              estimator_label, method_label, n_cores = detectCores() - 1) {
  
  # wrapper function for one simulation
  single_sim <- function(i) {
    df <- gen_data(n=1000, r = 1/4)
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

## Calculating the "true dataframes" 
true_CG_df <- data.frame(threshold = NB_true_quarter$threshold, NB_true = NB_true_quarter$NB_CG)
true_ATT_df <- data.frame(threshold = NB_true_quarter$threshold, NB_true = NB_true_quarter$NB_ATT)

## CG Coverage probability and CI Length #

CG_delta <- evaluate_coverage(
  n_sim = 2000,
  t.grd = t.grd,
  true_vals = true_CG_df,
  ci_function = function(y, a, pb, t.grd) {
    out <- delta_SE_CG(y, a, pb, t.grd)
    list(ci.lower = out$ci.lower, ci.upper = out$ci.upper)
  },
  estimator_label = "CG",
  method_label = "Delta"
)

View(CG_delta)

library(tictoc)

tic("Bayesian CG evaluation")

CG_bayes <- evaluate_coverage(
  n_sim = 2000,
  t.grd = t.grd,
  true_vals = true_CG_df,
  ci_function = function(y, a, pb, t.grd) {
    out <- bayesian_CG_Shared(y, a, pb, t.grd)
    list(ci.lower = out$ci.lower, ci.upper = out$ci.upper)
  },
  estimator_label = "CG",
  method_label = "Bayesian"
)

toc()

tic("Bootstrap CG evaluation")
 
CG_bootstrap <- evaluate_coverage(
  n_sim = 2000,
  t.grd = t.grd,
  true_vals = true_CG_df,
  ci_function = function(y, a, pb, t.grd) {
    out <- bootstrap_SE_CG(y, a, pb, t.grd, B = 1000)
    list(ci.lower = out$ci.lower, ci.upper = out$ci.upper)
  },
  estimator_label = "CG",
  method_label = "Bootstrap"
)

toc()

CG_combined <- bind_rows(CG_delta, CG_bayes, CG_bootstrap)


## ATT Coverage probability and CI Length #

ATT_delta <- evaluate_coverage(
  n_sim = 2000,
  t.grd = t.grd,
  true_vals = true_ATT_df,
  ci_function = function(y, a, pb, t.grd) {
    out <- delta_SE_ATT(y, a, pb, t.grd)
    list(ci.lower = out$ci.lower, ci.upper = out$ci.upper)
  },
  estimator_label = "ATT",
  method_label = "Delta"
)


tic("Bayes ATT evaluation")

ATT_bayes <- evaluate_coverage(
  n_sim = 2000,
  t.grd = t.grd,
  true_vals = true_ATT_df,
  ci_function = function(y, a, pb, t.grd) {
    out <- bayesian_ATT_Shared(y, a, pb, t.grd)
    list(ci.lower = out$nb.lower, ci.upper = out$nb.upper) 
  },
  estimator_label = "ATT",
  method_label = "Bayesian"
)

toc()

tic("Bootstrap ATT evaluation")

ATT_bootstrap <- evaluate_coverage(
  n_sim = 2000,
  t.grd = t.grd,
  true_vals = true_ATT_df,
  ci_function = function(y, a, pb, t.grd) {
    out <- bootstrap_SE_ATT(y, a, pb, t.grd, B = 1000)
    list(ci.lower = out$ci.lower, ci.upper = out$ci.upper)
  },
  estimator_label = "ATT",
  method_label = "Bootstrap"
)

toc()

ATT_combined <- bind_rows(ATT_delta, ATT_bayes, ATT_bootstrap)

View(ATT_combined)


