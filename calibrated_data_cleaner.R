######################################################################
# function: delta_sequence
# takes in:  length of the sequence, I
# returns: vector of length I equally spaced on the log scale
# note: for the ipsi case
delta_sequence_ipsi <- function(I) {
  log_lower <- -2.3
  log_upper <- 2.3
  
  # equally spaced on log scale
  log_delta_seq <- seq(log_lower, log_upper, length.out = I)  
  
  # back-transform to original scale
  delta_seq <- exp(log_delta_seq)  
  
  return(delta_seq)
}




######################################################################
# function: delta_sequence_static
# takes in:  length of the sequence, I
# returns: vector of length I equally spaced on the log scale
# note: for static case we have two deltas, see sim details. 
# This is lowkey cursed but it works lol. 
delta_sequence_static <- function(I) {
  delta_1 = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95)
  delta_2 = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95)
  
  delta_seq <- numeric(I)
  index <- 1
  
  for (i in 1:10) {
    for (j in 1:10) {
      delta_seq[index] = delta_1[i] + delta_2[j]
      index <- index + 1
    }
  }
  return(delta_seq)
}




##################################################################################
# function: p_X 
# takes in: covariates X, delta_seq
# returns matrix p_X containing all the delta values
p_X <- function(X, delta_seq) {
  I <- length(delta_seq)
  n <- nrow(X)
  p_X <- matrix(NA, nrow = n, ncol = I)
  
  for (i in 1:I) {
    d2_index <- i %% 10 
    if(d2_index == 0) {
      delta_2 <- 0.95
    } else{ 
      delta_2 <- -0.05 + 0.1 * d2_index
    }
    
    delta_1 <- delta_seq[i] - delta_2  
    
    for (j in 1:n) {
      if (X[j,1] <= 1) {
        p_X[j, i] = delta_1
      } else {
        p_X[j, i] = delta_2
      }
    }
    
  }
  return(p_X)
}







#######################################################################
# function: psi_true_ipsi
# takes in: a delta sequence value 
# returns: true incremental effect true_psi
# note: for a given delta finds the true expectation by using very large n
#       uses seed 42
psi_true_ipsi <- function(Delta_seq) {
  set.seed(42)
  I <- length(Delta_seq)
  true_psi <- numeric(I)
  
  n_true <- 1e6
  X_big <- matrix(rnorm(n_true * 4), ncol = 4)
  
  lin_pred <- -X_big[,1] + 0.5*X_big[,2] - 0.25*X_big[,3] - 0.1*X_big[,4]
  pi_true_big <- plogis(lin_pred)
  
  mu0 <- 200
  mu1 <- mu0 + 10 + 13.7 * (2*X_big[,1] + X_big[,2] + X_big[,3] + X_big[,4])
  
  # to itearate over delta sequence
  for (i in 1:I) {
    
    numerator <- Delta_seq[i] * pi_true_big * mu1 + (1 - pi_true_big)* mu0
    denominator <- Delta_seq[i]* pi_true_big + (1 - pi_true_big)
    
    true_psi[i] <- mean(numerator / denominator)
    
  }
  return(true_psi)
}




#######################################################################
# function: psi_true_static
# takes in: a delta sequence value 
# returns: true incremental effect true_psi
# note: for a given delta finds the true expectation by using very large n
#       uses seed 42
psi_true_static <- function(delta_seq) {
  set.seed(42)
  I <- length(delta_seq)
  true_psi <- numeric(I)
  
  n_true <- 1e6
  X_big <- matrix(rnorm(n_true * 4), ncol = 4)
  
  lin_pred <- -X_big[,1] + 0.5*X_big[,2] - 0.25*X_big[,3] - 0.1*X_big[,4]
  pi_true_big <- plogis(lin_pred)
  
  mu0 <- 200
  mu1 <- mu0 + 10 + 13.7 * (2*X_big[,1] + X_big[,2] + X_big[,3] + X_big[,4])
  p_mat = p_X(X_big, delta_seq)
  
  # to itearate over delta sequence
  for (i in 1:I) {
    true_psi[i] <- mean(mu1 * p_mat[,i] + mu0 * (1-p_mat[,i]))
  }
  return(true_psi)
}




#########################################################################
# function: bias_est
# takes in: est_matrix, true_vec, I, J
# returns: bias
bias_est <- function(est_matrix, true_vec, I, J) {
  
  # holds the absolute value of internal loop for averaging
  ans_vec = numeric(I)
  
  
  
  for (i in 1:I) {
    est_vec <- est_matrix[, i]
    diff_vec = est_vec - true_vec[i]
    sum = 0
    for (j in 1:J) {
      sum = sum + diff_vec[j]
    }
    ans_vec[i] <- abs(sum / J)
  }
  bias = mean(ans_vec)
  return(bias)
}




##########################################################################
# function: RMSE_est
# takes in: est_vec, true_vec, I, J, n
# returns: RMSE
RMSE_est <- function(est_matrix, true_vec, I, J, n) {
  
  # holds the absolute value of internal loop for averaging
  ans_vec = numeric(I)
  
  
  
  
  for (i in 1:I) {
    est_vec <- est_matrix[, i]
    diff_vec = est_vec - true_vec[i]
    sum = 0
    for (j in 1:J) {
      sum = sum + diff_vec[j]^2
    }
    ans_vec[i] <- sqrt(sum / J)
  }
  RMSE = sqrt(n) * mean(ans_vec)
  return(RMSE)
}

################### extract sims
data_cleaner <- function(n, I, J, Version, Size) {
  ipsi_delta_sequence <- delta_sequence_ipsi(I=I)
  ipsi_true_psi <- psi_true_ipsi(Delta_seq = ipsi_delta_sequence)
  
  static_delta_sequence <- delta_sequence_static(I=I)
  static_true_psi <- psi_true_static(delta_seq = static_delta_sequence)
  
  all_sims <- vector("list", J)  # preallocate a list
  
  for (j in 1:J) {
    filename <- paste0(Version, Size, "_sim", j, ".rds")
    all_sims[[j]] <- readRDS(filename)
  }
  
  # data storage setup
  
  # ipsi reg
  ipsi_reg_psi_hat <- matrix(NA, nrow = J, ncol = I)
  ipsi_reg_efficient_psi_hat <- matrix(NA, nrow = J, ncol = I)
  
  ipsi_reg_uniform_coverage <- numeric(J)
  ipsi_reg_pointwise_coverage <- matrix(NA, nrow = J, ncol = I)
  ipsi_reg_uniform_length <- matrix(NA, nrow = J, ncol = I)
  ipsi_reg_pointwise_length <- matrix(NA, nrow = J, ncol = I)
  
  ipsi_reg_efficient_uniform_coverage <- numeric(J)
  ipsi_reg_efficient_pointwise_coverage <- matrix(NA, nrow = J, ncol = I)
  ipsi_reg_efficient_uniform_length <- matrix(NA, nrow = J, ncol = I)
  ipsi_reg_efficient_pointwise_length <- matrix(NA, nrow = J, ncol = I)
  
  # ipsi trans
  ipsi_trans_psi_hat <-  matrix(NA, nrow = J, ncol = I)
  ipsi_trans_efficient_psi_hat <- matrix(NA, nrow = J, ncol = I)
  
  ipsi_trans_uniform_coverage <- numeric(J)
  ipsi_trans_pointwise_coverage <- matrix(NA, nrow = J, ncol = I)
  ipsi_trans_uniform_length <- matrix(NA, nrow = J, ncol = I)
  ipsi_trans_pointwise_length <- matrix(NA, nrow = J, ncol = I)
  
  ipsi_trans_efficient_uniform_coverage <- numeric(J)
  ipsi_trans_efficient_pointwise_coverage <- matrix(NA, nrow = J, ncol = I)
  ipsi_trans_efficient_uniform_length <- matrix(NA, nrow = J, ncol = I)
  ipsi_trans_efficient_pointwise_length <- matrix(NA, nrow = J, ncol = I)
  
  # static reg
  static_reg_psi_hat <- matrix(NA, nrow = J, ncol = I)
  static_reg_efficient_psi_hat <- matrix(NA, nrow = J, ncol = I)
  
  static_reg_uniform_coverage <- numeric(J)
  static_reg_pointwise_coverage <- matrix(NA, nrow = J, ncol = I)
  static_reg_uniform_length <- matrix(NA, nrow = J, ncol = I)
  static_reg_pointwise_length <- matrix(NA, nrow = J, ncol = I)
  
  static_reg_efficient_uniform_coverage <- numeric(J)
  static_reg_efficient_pointwise_coverage <- matrix(NA, nrow = J, ncol = I)
  static_reg_efficient_uniform_length <- matrix(NA, nrow = J, ncol = I)
  static_reg_efficient_pointwise_length <- matrix(NA, nrow = J, ncol = I)
  
  # static trans
  static_trans_psi_hat <-  matrix(NA, nrow = J, ncol = I)
  static_trans_efficient_psi_hat <- matrix(NA, nrow = J, ncol = I)
  
  static_trans_uniform_coverage <- numeric(J)
  static_trans_pointwise_coverage <- matrix(NA, nrow = J, ncol = I)
  static_trans_uniform_length <- matrix(NA, nrow = J, ncol = I)
  static_trans_pointwise_length <- matrix(NA, nrow = J, ncol = I)
  
  static_trans_efficient_uniform_coverage <- numeric(J)
  static_trans_efficient_pointwise_coverage <- matrix(NA, nrow = J, ncol = I)
  static_trans_efficient_uniform_length <- matrix(NA, nrow = J, ncol = I)
  static_trans_efficient_pointwise_length <- matrix(NA, nrow = J, ncol = I)
  
  # extraction
  for (j in 1:J) {
    # ipsi reg
    ipsi_reg_psi_hat[j,] <- all_sims[[j]]$ipsi_reg_psi_hat
    ipsi_reg_efficient_psi_hat[j,] <- all_sims[[j]]$ipsi_reg_efficient_psi_hat
    
    ipsi_reg_uniform_coverage[j] <- all_sims[[j]]$ipsi_reg_coverage$uniform_coverage
    ipsi_reg_pointwise_coverage[j, ] <- all_sims[[j]]$ipsi_reg_coverage$pointwise_coverage
    ipsi_reg_uniform_length[j,] <- all_sims[[j]]$ipsi_reg_coverage$uniform_length
    ipsi_reg_pointwise_length[j, ] <- all_sims[[j]]$ipsi_reg_coverage$pointwise_length
    
    ipsi_reg_efficient_uniform_coverage[j] <- all_sims[[j]]$ipsi_reg_efficient_coverage$uniform_coverage
    ipsi_reg_efficient_pointwise_coverage[j, ] <- all_sims[[j]]$ipsi_reg_efficient_coverage$pointwise_coverage
    ipsi_reg_efficient_uniform_length[j,] <- all_sims[[j]]$ipsi_reg_efficient_coverage$uniform_length
    ipsi_reg_efficient_pointwise_length[j, ] <- all_sims[[j]]$ipsi_reg_efficient_coverage$pointwise_length
    
    # ipsi trans
    ipsi_trans_psi_hat[j,] <- all_sims[[j]]$ipsi_trans_psi_hat
    ipsi_trans_efficient_psi_hat[j,] <- all_sims[[j]]$ipsi_trans_efficient_psi_hat
    
    ipsi_trans_uniform_coverage[j] <- all_sims[[j]]$ipsi_trans_coverage$uniform_coverage
    ipsi_trans_pointwise_coverage[j, ] <- all_sims[[j]]$ipsi_trans_coverage$pointwise_coverage
    ipsi_trans_uniform_length[j, ] <- all_sims[[j]]$ipsi_trans_coverage$uniform_length
    ipsi_trans_pointwise_length[j, ] <- all_sims[[j]]$ipsi_trans_coverage$pointwise_length
    
    ipsi_trans_efficient_uniform_coverage[j] <- all_sims[[j]]$ipsi_trans_efficient_coverage$uniform_coverage
    ipsi_trans_efficient_pointwise_coverage[j, ] <- all_sims[[j]]$ipsi_trans_efficient_coverage$pointwise_coverage
    ipsi_trans_efficient_uniform_length[j, ] <- all_sims[[j]]$ipsi_trans_efficient_coverage$uniform_length
    ipsi_trans_efficient_pointwise_length[j, ] <- all_sims[[j]]$ipsi_trans_efficient_coverage$pointwise_length
    
    # static reg
    static_reg_psi_hat[j,] <- all_sims[[j]]$static_reg_psi_hat
    static_reg_efficient_psi_hat[j,] <- all_sims[[j]]$static_reg_efficient_psi_hat
    
    static_reg_uniform_coverage[j] <- all_sims[[j]]$static_reg_coverage$uniform_coverage
    static_reg_pointwise_coverage[j, ] <- all_sims[[j]]$static_reg_coverage$pointwise_coverage
    static_reg_uniform_length[j,] <- all_sims[[j]]$static_reg_coverage$uniform_length
    static_reg_pointwise_length[j, ] <- all_sims[[j]]$static_reg_coverage$pointwise_length
    
    static_reg_efficient_uniform_coverage[j] <- all_sims[[j]]$static_reg_efficient_coverage$uniform_coverage
    static_reg_efficient_pointwise_coverage[j, ] <- all_sims[[j]]$static_reg_efficient_coverage$pointwise_coverage
    static_reg_efficient_uniform_length[j,] <- all_sims[[j]]$static_reg_efficient_coverage$uniform_length
    static_reg_efficient_pointwise_length[j, ] <- all_sims[[j]]$static_reg_efficient_coverage$pointwise_length
    
    # static trans
    static_trans_psi_hat[j,] <- all_sims[[j]]$static_trans_psi_hat
    static_trans_efficient_psi_hat[j,] <- all_sims[[j]]$static_trans_efficient_psi_hat
    
    static_trans_uniform_coverage[j] <- all_sims[[j]]$static_trans_coverage$uniform_coverage
    static_trans_pointwise_coverage[j, ] <- all_sims[[j]]$static_trans_coverage$pointwise_coverage
    static_trans_uniform_length[j, ] <- all_sims[[j]]$static_trans_coverage$uniform_length
    static_trans_pointwise_length[j, ] <- all_sims[[j]]$static_trans_coverage$pointwise_length
    
    static_trans_efficient_uniform_coverage[j] <- all_sims[[j]]$static_trans_efficient_coverage$uniform_coverage
    static_trans_efficient_pointwise_coverage[j, ] <- all_sims[[j]]$static_trans_efficient_coverage$pointwise_coverage
    static_trans_efficient_uniform_length[j, ] <- all_sims[[j]]$static_trans_efficient_coverage$uniform_length
    static_trans_efficient_pointwise_length[j, ] <- all_sims[[j]]$static_trans_efficient_coverage$pointwise_length
  }
  
  # ipsi
  # bias
  ipsi_reg_bias <- bias_est(est_matrix = ipsi_reg_psi_hat, true_vec = ipsi_true_psi, I= I, J=J)
  ipsi_reg_efficient_bias <- bias_est(est_matrix = ipsi_reg_efficient_psi_hat, true_vec = ipsi_true_psi, I= I, J=J)
  
  # rmse
  ipsi_reg_rmse <- RMSE_est(est_matrix = ipsi_reg_psi_hat, true_vec = ipsi_true_psi, I= I, J=J, n=n)
  ipsi_reg_efficient_rmse <- RMSE_est(est_matrix = ipsi_reg_efficient_psi_hat, true_vec = ipsi_true_psi, I= I, J=J, n=n)
  
  # trans bias
  ipsi_trans_bias <- bias_est(est_matrix = ipsi_trans_psi_hat, true_vec = ipsi_true_psi, I= I, J=J)
  ipsi_trans_efficient_bias <- bias_est(est_matrix = ipsi_trans_efficient_psi_hat, true_vec = ipsi_true_psi, I= I, J=J)
  
  # trans rmse 
  ipsi_trans_rmse <- RMSE_est(est_matrix = ipsi_trans_psi_hat, true_vec = ipsi_true_psi, I= I, J=J, n=n)
  ipsi_trans_efficient_rmse <- RMSE_est(est_matrix = ipsi_trans_efficient_psi_hat, true_vec = ipsi_true_psi, I= I, J=J, n=n)
  
  # static
  # bias
  static_reg_bias <- bias_est(est_matrix = static_reg_psi_hat, true_vec = static_true_psi, I= I, J=J)
  static_reg_efficient_bias <- bias_est(est_matrix = static_reg_efficient_psi_hat, true_vec = static_true_psi, I= I, J=J)
  
  # rmse
  static_reg_rmse <- RMSE_est(est_matrix = static_reg_psi_hat, true_vec = static_true_psi, I= I, J=J, n=n)
  static_reg_efficient_rmse <- RMSE_est(est_matrix = static_reg_efficient_psi_hat, true_vec = static_true_psi, I= I, J=J, n=n)
  
  # trans bias
  static_trans_bias <- bias_est(est_matrix = static_trans_psi_hat, true_vec = static_true_psi, I= I, J=J)
  static_trans_efficient_bias <- bias_est(est_matrix = static_trans_efficient_psi_hat, true_vec = static_true_psi, I= I, J=J)
  
  # trans rmse 
  static_trans_rmse <- RMSE_est(est_matrix = static_trans_psi_hat, true_vec = static_true_psi, I= I, J=J, n=n)
  static_trans_efficient_rmse <- RMSE_est(est_matrix = static_trans_efficient_psi_hat, true_vec = static_true_psi, I= I, J=J, n=n)
  
  return(list(
    ipsi_reg_uniform_coverage = ipsi_reg_uniform_coverage,
    ipsi_reg_pointwise_coverage = ipsi_reg_pointwise_coverage,
    ipsi_reg_uniform_length = ipsi_reg_uniform_length,
    ipsi_reg_pointwise_length = ipsi_reg_pointwise_length,
    
    ipsi_reg_efficient_uniform_coverage = ipsi_reg_efficient_uniform_coverage,
    ipsi_reg_efficient_pointwise_coverage = ipsi_reg_efficient_pointwise_coverage,
    ipsi_reg_efficient_uniform_length = ipsi_reg_efficient_uniform_length,
    ipsi_reg_efficient_pointwise_length = ipsi_reg_efficient_pointwise_length,
    
    ipsi_reg_bias = ipsi_reg_bias, 
    ipsi_reg_efficient_bias = ipsi_reg_efficient_bias, 
    ipsi_reg_rmse = ipsi_reg_rmse, 
    ipsi_reg_efficient_rmse = ipsi_reg_efficient_rmse,
    
    ipsi_trans_uniform_coverage = ipsi_trans_uniform_coverage,
    ipsi_trans_pointwise_coverage = ipsi_trans_pointwise_coverage,
    ipsi_trans_uniform_length = ipsi_trans_uniform_length,
    ipsi_trans_pointwise_length = ipsi_trans_pointwise_length,
    
    ipsi_trans_efficient_uniform_coverage = ipsi_trans_efficient_uniform_coverage,
    ipsi_trans_efficient_pointwise_coverage = ipsi_trans_efficient_pointwise_coverage,
    ipsi_trans_efficient_uniform_length = ipsi_trans_efficient_uniform_length,
    ipsi_trans_efficient_pointwise_length = ipsi_trans_efficient_pointwise_length,
    
    ipsi_trans_bias = ipsi_trans_bias, 
    ipsi_trans_efficient_bias = ipsi_trans_efficient_bias,
    ipsi_trans_rmse = ipsi_trans_rmse,
    ipsi_trans_efficient_rmse = ipsi_trans_efficient_rmse,
    
    static_reg_uniform_coverage = static_reg_uniform_coverage,
    static_reg_pointwise_coverage = static_reg_pointwise_coverage,
    static_reg_uniform_length = static_reg_uniform_length,
    static_reg_pointwise_length = static_reg_pointwise_length,
    
    static_reg_efficient_uniform_coverage = static_reg_efficient_uniform_coverage,
    static_reg_efficient_pointwise_coverage = static_reg_efficient_pointwise_coverage,
    static_reg_efficient_uniform_length = static_reg_efficient_uniform_length,
    static_reg_efficient_pointwise_length = static_reg_efficient_pointwise_length,
    
    static_reg_bias = static_reg_bias, 
    static_reg_efficient_bias = static_reg_efficient_bias, 
    static_reg_rmse = static_reg_rmse, 
    static_reg_efficient_rmse = static_reg_efficient_rmse,
    
    static_trans_uniform_coverage = static_trans_uniform_coverage,
    static_trans_pointwise_coverage = static_trans_pointwise_coverage,
    static_trans_uniform_length = static_trans_uniform_length,
    static_trans_pointwise_length = static_trans_pointwise_length,
    
    static_trans_efficient_uniform_coverage =  static_trans_efficient_uniform_coverage,
    static_trans_efficient_pointwise_coverage = static_trans_efficient_pointwise_coverage,
    static_trans_efficient_uniform_length = static_trans_efficient_uniform_length,
    static_trans_efficient_pointwise_length = static_trans_efficient_pointwise_length,
    
    static_trans_bias = static_trans_bias, 
    static_trans_efficient_bias = static_trans_efficient_bias,
    static_trans_rmse = static_trans_rmse,
    static_trans_efficient_rmse = static_trans_efficient_rmse
  ))
}
Version = "softbcf"
Size = "n3"

data <- data_cleaner(n = 5000, I = 100, J = 500, Version, Size)

filename <- paste0(Version, "_",Size, "_completed_data.rds")
saveRDS(data, file=filename)

