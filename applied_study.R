############################################################
# Load packages
############################################################

library(nhanesA)
library(dplyr)
library(tidyr)
library(SoftBart)

###########################################################################################
# function: normalize01
# note: quantile normalization to [0,1]
normalize01 <- function(x) {
  # scales by x and evaluates at x
  ecdf(x)(x)   # empirical CDF transform
}




#######################################################################
# function: fit_probit
# takes in: covariates X, outcome Y, num_tree, num_iter
# returns: (num_iter x n) posterior draws matrix pi_train, 
#          burn-in removed pi_hat.
# note: the burn-in calculation implicitly assumes n is even
fit_probit <- function(X, Y, num_tree, num_iter=4000) {
  hypers <- Hypers(X, Y, k = 1/6, num_tree = num_tree, sigma_hat = 1)
  opts   <- Opts(update_sigma = FALSE)
  probit_forest <- MakeForest(hypers, opts)
  
  # store posterior draws of pi(x) = Pr(Y=1|X)
  pi_train <- matrix(nrow = num_iter, ncol = nrow(X))
  
  # initialize
  r <- probit_forest$do_predict(X)
  upper <- ifelse(Y == 0, 0, Inf)
  lower <- ifelse(Y == 0, -Inf, 0)
  Z <- truncnorm::rtruncnorm(n = length(Y), a = lower, b = upper, mean = r, sd = 1)
  
  for(i in 1:num_iter) {
    r <- probit_forest$do_gibbs(X, Z, X, 1)
    Z <- truncnorm::rtruncnorm(n = length(Y), a = lower, b = upper, mean = r, sd = 1)
    pi_train[i, ] <- pnorm(r)  # posterior draw of nuisance function
  }
  pi_hat <- colMeans(pi_train[(num_iter/2 +1): num_iter, ])
  return(list(pi_train = pi_train, pi_hat = pi_hat))
}




#####################################################################################################
# function: fit_softbart
# takes in: outcome y, treatment Z, covariates X, num_iter, num_burn, num_tree, verbose (for debugging)
# returns: 
fit_softbart <- function(y, Z, X, num_iter = 4000, num_burn = 0, num_tree = 200, verbose = TRUE) {
  n <- length(y)
  
  # Combine treatment and covariates
  X_full <- cbind(Z = Z, X)
  
  # Potential outcome covariate matrices
  X0 <- cbind(Z = rep(0, n), X)
  X1 <- cbind(Z = rep(1, n), X)
  
  # Set up hyperparameters and MCMC options
  hypers <- Hypers(X = X_full, Y = y, num_tree = num_tree)
  opts <- Opts(num_burn = num_burn, num_save = num_iter)
  
  # Run SoftBART
  fit <- softbart(
    X = X_full,
    Y = y,
    X_test = rbind(X0, X1),
    hypers = hypers,
    opts = opts,
    verbose = verbose
  )
  
  # Extract posterior predictions
  mu_test <- fit$y_hat_test
  mu_train <- fit$y_hat_train
  
  # Sanity check: remove missing or NaN rows
  good_draws <- which(rowSums(is.na(mu_test)) == 0)
  if (length(good_draws) < nrow(mu_test)) {
    warning(sprintf("Removed %d draws with NAs", nrow(mu_test) - length(good_draws)))
  }
  mu_test <- mu_test[good_draws, , drop = FALSE]
  mu_train <- mu_train[good_draws, , drop = FALSE]
  sigma <- fit$sigma[good_draws]
  
  # Split posterior predictions into mu0 and mu1 parts
  mu0_post <- mu_test[, 1:n, drop = FALSE]
  mu1_post <- mu_test[, (n + 1):(2 * n), drop = FALSE]
  
  return(list(
    mu0_hat = mu0_post,
    mu1_hat = mu1_post,
    f = mu_train,
    sigma = sigma
  ))
}


###################################################################################
# function: bayes_boot
# takes in: covar X, treatment Z, outcome y, number of bootstrap iterations B,
#           list from nuisance fit function, delta sequence, and intervention (either "ipsi" or "static")
# returns: (B x I) matrices psi_delta_matrix and efficient_influence_matrix 
#.          for the plug-in and one-step estimators respectively.
bayes_boot <- function(X, Z, y, B, nuisance_fit, delta_seq) {
  
  n <- length(Z) 
  I <- length(delta_seq) 
  m_t_matrix <- matrix(NA, nrow = B, ncol = I)
  psi_delta_matrix <- matrix(NA, nrow = B, ncol = I)
  efficient_influence_matrix <- matrix(NA, nrow = B, ncol = I)
  
  
  
  pi_hat <- nuisance_fit$pi_hat_distribution
  mu0_hat <- nuisance_fit$mu0_hat_distribution
  mu1_hat <- nuisance_fit$mu1_hat_distribution
  
  
  # burn point
  num_iter <- 4000
  burn <- floor(num_iter / 2)
  size <- num_iter - burn
  
    for (b in 1:B) {
      # dirchlet distribution is the same as the standardized exponential distribution
      
      
      random_exp <- rexp(n)
      
      random_dirch <- random_exp / sum(random_exp)
      
      pi_hat_sample <- sample(1:size, 1)
      mu_sample <- sample(1:size, 1)
      
      # grabs a random pi_hat, mu1_hat and mu0_hat to account for variability in their models
      pi_hat_draw <- pi_hat[pi_hat_sample, ]
      mu0_hat_draw <- mu0_hat[mu_sample, ]
      mu1_hat_draw <- mu1_hat[mu_sample, ]
      
      
      
      for (i in 1:I) {
        delta <- delta_seq[i]
        # denominator
        denom <- delta * pi_hat_draw + (1 - pi_hat_draw)
        ratio <- (delta * pi_hat_draw * mu1_hat_draw + (1 - pi_hat_draw) * mu0_hat_draw) / denom
        
        psi_delta_matrix[b, i] <- sum(ratio * random_dirch)
        
        efficient_influence_value <- ipsi_efficient_influence(Z=Z, y=y, pi_hat=pi_hat_draw,
                                                              delta=delta, m_t_1=mu1_hat_draw,m_t_0=mu0_hat_draw, ratio = ratio)
        
        efficient_influence_matrix[b, i] <- mean(ratio) + sum(efficient_influence_value * random_dirch)                                                  
      }
      
      
      ## progress update every 100 iterations
      if (b %% 100 == 0) {
        cat("Finished bootstrap", b, "of", B, "\n")
        flush.console() 
      }
      
    }
  
  
  
  return(list(psi_delta_matrix = psi_delta_matrix, efficient_influence_matrix = efficient_influence_matrix))
}




##################################################################################
# function: ipsi_efficient_influence 
# takes in: Z, y, pi_hat_draw, delta, m_t, t
# returns efficient influence function estimate psi_delta
# note: this algebraic form is from the longitudinal form definition simplified 
#       to T = 1.
ipsi_efficient_influence <- function(Z, y, pi_hat, delta, m_t_1, m_t_0, ratio) {
  score_num <- Z *(1- pi_hat) - (1-Z) * delta * pi_hat
  score_denom <- delta / (1-delta)
  score_term <- score_num / score_denom
  
  shift_num <- delta * pi_hat * m_t_1 + (1-pi_hat) * m_t_0
  shift_denom <- delta * pi_hat + 1 - pi_hat
  shift_term <- shift_num / shift_denom
  
  cum_weight_num <- delta * Z + 1 - Z 
  cum_weight_denom <- delta * pi_hat + 1 - pi_hat
  cum_weight_term <- cum_weight_num/ cum_weight_denom
  
  y_num <- (delta * Z + 1 -Z) * y
  y_denom <- delta *pi_hat + 1 - pi_hat
  y_term <- y_num / y_denom
  
  phi_hat <- score_term * shift_term * cum_weight_term + y_term - mean(ratio)
  return(phi_hat)
}

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
# function: coverage
# takes in: (B x I) 
# returns if 95% credible interval contains true psi
bayesian_coverage <- function(est_matrix) {
  psi_hat <- colMeans(est_matrix)  # posterior mean at each delta
  
  B <- nrow(est_matrix)
  I <- ncol(est_matrix)
  
  # Posterior standard deviation at each delta (manual)
  sigma_hat <- numeric(I)
  for (i in 1:I) {
    sigma_hat[i] <- sqrt(sum((est_matrix[, i] - psi_hat[i])^2) / (B - 1))
  }
  
  max_vals <- numeric(B)
  
  # Pointwise coverage and length
  eff_ll <- numeric(I)
  eff_ul <- numeric(I)
  
  
  for (i in 1:I) {
    eff_ll[i] <- quantile(est_matrix[, i], 0.025)
    eff_ul[i] <- quantile(est_matrix[, i], 0.975)
  }
  
  # Uniform coverage (studentized)
  for (b in 1:B) {
    max_vals[b] <- max(abs(est_matrix[b, ] - psi_hat) / sigma_hat)
  }
  
  crit_val <- quantile(max_vals, 0.95)
  
  uniform_ll <- psi_hat - crit_val * sigma_hat
  uniform_ul <- psi_hat + crit_val * sigma_hat
  uniform_length <- uniform_ul - uniform_ll
  
  
  
  return(list(
    pointwise_lower = eff_ll,
    pointwise_upper = eff_ul,
    uniform_lower = uniform_ll, 
    uniform_upper = uniform_ul
  ))
}



###################################################################################
# ********************************************************************************#
###################################################################################
# End of Functions





############################################################
# Download NHANES 2017–2018 data
############################################################

# Demographics
demo <- nhanes("DEMO_J")

# Prescription medications
rx <- nhanes("RXQ_RX_J")

# Body measures (BMI)
bmx <- nhanes("BMX_J")

# LDL cholesterol
trig <- nhanes("TRIGLY_J")   # contains LBDLDL

# Diabetes questionnaire
diq <- nhanes("DIQ_J")

# Smoking questionnaire
smq <- nhanes("SMQ_J")

# Blood pressure questionnaire (hypertension)
bpq <- nhanes("BPQ_J")

############################################################
# Define statin use indicator
############################################################

statins <- c(
  "ATORVASTATIN",
  "SIMVASTATIN",
  "ROSUVASTATIN",
  "PRAVASTATIN",
  "LOVASTATIN",
  "FLUVASTATIN",
  "PITAVASTATIN"
)

statin_use <- rx %>%
  mutate(
    drug_upper = toupper(RXDDRUG),
    statin = drug_upper %in% statins
  ) %>%
  group_by(SEQN) %>%
  summarise(
    statin = as.integer(any(statin)),
    .groups = "drop"
  )

############################################################
# Merge all datasets
############################################################

data <- demo %>%
  left_join(bmx, by = "SEQN") %>%
  left_join(trig, by = "SEQN") %>%
  left_join(diq, by = "SEQN") %>%
  left_join(smq, by = "SEQN") %>%
  left_join(bpq, by = "SEQN") %>%
  left_join(statin_use, by = "SEQN")

############################################################
# Create treatment, outcome, and covariates
############################################################



############################################################
# Preprocessing for modeling
############################################################

# 1. Filter adults and remove rows with essential missing values
data_clean <- data %>%
  filter(RIDAGEYR >= 18) %>%
  mutate(
    # Numeric covariates
    age = RIDAGEYR,
    bmi = BMXBMI,
    
    # Binary covariates
    sex = ifelse(RIAGENDR == "Male" | RIAGENDR == 1, 1, 0),
    diabetes = ifelse(DIQ010 == "Yes", 1,
                      ifelse(DIQ010 == "No", 0, NA)),
    smoker = ifelse(SMQ020 == "Yes", 1,
                    ifelse(SMQ020 == "No", 0, NA)),
    hypertension = ifelse(BPQ020 == "Yes", 1,
                          ifelse(BPQ020 == "No", 0, NA)),
    
    # Categorical covariate
    race = factor(RIDRETH3),
    
    # Outcome and treatment
    ldl = LBDLDL
    # statin is already 0/1 from statin_use
  ) %>%
  select(age, sex, bmi, diabetes, smoker, hypertension, race, ldl, statin) %>%
  drop_na(age, sex, bmi, diabetes, smoker, hypertension, race, ldl, statin)


############################################################
# 2. Treatment and outcome vectors
############################################################

A <- as.vector(data_clean$statin)
y <- as.vector(data_clean$ldl)


############################################################
# 3. Normalize numeric covariates (age, bmi)
############################################################

normalize01 <- function(x) ecdf(x)(x)

X_num <- data_clean %>%
  select(age, bmi, sex, diabetes, smoker, hypertension) %>%
  mutate(
    age = normalize01(age),
    bmi = normalize01(bmi)
  ) %>%
  as.matrix()


############################################################
# 4. Convert categorical variable race into dummy variables
############################################################

X_cat <- model.matrix(~ race - 1, data = data_clean)  # one-hot encoding

############################################################
# 5. Combine numeric + categorical covariates
############################################################

X_norm <- cbind(X_num, X_cat)


# calling functions to fit the models
# Fit probit model for statin use
probit_fit <- fit_probit(X = X_norm, Y = A, num_tree = 50, num_iter = 4000)

pi_hat_distribution <- probit_fit$pi_train
pi_hat <- probit_fit$pi_hat




# Fit SoftBart outcome regression with treatment included
outcome_fit <- fit_softbart(y = y, Z = A, X = X_norm, num_iter = 4000, num_burn = 0, num_tree = 200)

# Extract posterior draws and potential outcomes
mu0_hat_burned <- outcome_fit$mu0_hat[2001:4000, ]
mu1_hat_burned <- outcome_fit$mu1_hat[2001:4000, ]


# used for trace plots
mu0_hat <- outcome_fit$mu0_hat
mu1_hat <- outcome_fit$mu1_hat


nuisance_fit <- list(pi_hat = pi_hat, 
                     pi_hat_distribution = pi_hat_distribution, 
                     mu0_hat_distribution = mu0_hat_burned, 
                     mu1_hat_distribution = mu1_hat_burned)
I = 100

ipsi_delta_seq <-delta_sequence_ipsi((I=I))

ipsi_psi_matrix <- bayes_boot(X = X_norm, Z = A, y = y, B = 10000, nuisance_fit = nuisance_fit,
                                  delta_seq = ipsi_delta_seq)
ipsi_psi_hat <- colMeans(ipsi_psi_matrix$efficient_influence_matrix)

###############################################################################################



library(ggplot2)
library(reshape2)

# pick an individual
i <- 1  # ith person

# Extract posterior draws
pi_draws_i <- probit_fit$pi_train[, i]
mu1_draws_i <- mu1_hat[, i]
mu0_draws_i <- mu0_hat[, i]

# Combine into a data frame for plotting
mcmc_df <- data.frame(
  Iteration = 1:length(pi_draws_i),
  mu1_hat = mu1_draws_i,
  mu0_hat = mu0_draws_i
)

mcmc_df_long <- melt(mcmc_df, id.vars = "Iteration")

ggplot(mcmc_df_long, aes(x = Iteration, y = value, color = variable)) +
  geom_line(alpha = 0.7) +
  labs(title = paste("MCMC Trace for Individual", i),
       y = "Posterior Draw",
       color = "Parameter") +
  theme_minimal()

# Trace plot for pi only
i <- 1  # ith person
pi_draws_i <- probit_fit$pi_train[, i]

ggplot(data.frame(Iteration = 1:length(pi_draws_i), pi_hat = pi_draws_i),
       aes(x = Iteration, y = pi_hat)) +
  geom_line(color = "darkgreen") +
  labs(title = paste("MCMC Trace for pi_hat (Individual", i, ")"),
       x = "Iteration", y = expression(pi[hat])) +
  theme_minimal()

#####################################################################

sigma_draws <- outcome_fit$sigma

# Histogram + density
ggplot(data.frame(sigma = sigma_draws), aes(x = sigma)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", alpha = 0.6) +
  geom_density(color = "red", size = 1) +
  labs(title = "Posterior Distribution of Sigma", x = "Sigma", y = "Density") +
  theme_minimal()


# QQ-plot to check approximate normality
ggplot(data.frame(sigma = sigma_draws), aes(sample = sigma)) +
  stat_qq(color = "purple") +
  stat_qq_line(color = "red") +
  labs(title = "QQ-Plot of Posterior Sigma",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_minimal()

# Trace plot
ggplot(data.frame(Iteration = 1:length(sigma_draws), sigma = sigma_draws),
       aes(x = Iteration, y = sigma)) +
  geom_line(color = "blue") +
  labs(title = "Sigma Trace Plot", x = "Iteration", y = "Sigma") +
  theme_minimal()




pi_hat <- probit_fit$pi_hat  # posterior mean

ggplot(data.frame(pi_hat = pi_hat, A = factor(A)),
       aes(x = pi_hat, fill = A)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = "Propensity Score Overlap by Treatment",
       x = "pi_hat (Pr(A=1|X))", fill = "Treatment") +
  theme_minimal()

###############################################################################
# eif version
confidence <- bayesian_coverage(ipsi_psi_matrix$efficient_influence_matrix)


# Create summary dataframe for EIF estimator
delta_df <- data.frame(
  delta = ipsi_delta_seq,
  psi_mean = colMeans(ipsi_psi_matrix$efficient_influence_matrix),
  psi_lower = confidence$pointwise_lower,
  psi_upper = confidence$pointwise_upper,
  psi_lower_uniform = confidence$uniform_lower,
  psi_upper_uniform = confidence$uniform_upper
)


ribbons_df <- data.frame(
  delta = rep(delta_df$delta, 2),
  lower = c(delta_df$psi_lower, delta_df$psi_lower_uniform),
  upper = c(delta_df$psi_upper, delta_df$psi_upper_uniform),
  type = rep(c("Pointwise 95% CI", "Uniform 95% CI"), each = nrow(delta_df))
)

# Plot
ggplot() +
  geom_ribbon(data = ribbons_df, 
              aes(x = delta, ymin = lower, ymax = upper, fill = type),
              alpha = 0.3) +
  geom_line(data = delta_df, aes(x = delta, y = psi_mean, color = "Posterior Mean"), size = 1) +
  scale_fill_manual(values = c("Pointwise 95% CI" = "blue", "Uniform 95% CI" = "red")) +
  scale_color_manual(values = c("Posterior Mean" = "black")) +
  scale_x_log10() +
  labs(title = "Estimated Shift Effect vs Delta, EIF",
       x = "Delta (log scale)",
       y = expression(psi[hat]),
       fill = "Interval Type",
       color = "") +
  theme_minimal()



###############################################################################
# plug-in version
plug_confidence <- bayesian_coverage(ipsi_psi_matrix$psi_delta_matrix)


plug_delta_df <- data.frame(
  delta = ipsi_delta_seq,
  psi_mean = colMeans(ipsi_psi_matrix$psi_delta_matrix),
  psi_lower = plug_confidence$pointwise_lower,
  psi_upper = plug_confidence$pointwise_upper,
  psi_lower_uniform = plug_confidence$uniform_lower,
  psi_upper_uniform = plug_confidence$uniform_upper
)


plug_ribbons_df <- data.frame(
  delta = rep(plug_delta_df$delta, 2),
  lower = c(plug_delta_df$psi_lower, plug_delta_df$psi_lower_uniform),
  upper = c(plug_delta_df$psi_upper, plug_delta_df$psi_upper_uniform),
  type = rep(c("Pointwise 95% CI", "Uniform 95% CI"), each = nrow(plug_delta_df))
)

# Plot
ggplot() +
  geom_ribbon(data = plug_ribbons_df, 
              aes(x = delta, ymin = lower, ymax = upper, fill = type),
              alpha = 0.3) +
  geom_line(data = plug_delta_df, aes(x = delta, y = psi_mean, color = "Posterior Mean"), size = 1) +
  scale_fill_manual(values = c("Pointwise 95% CI" = "blue", "Uniform 95% CI" = "red")) +
  scale_color_manual(values = c("Posterior Mean" = "black")) +
  scale_x_log10() +
  labs(title = "Estimated Shift Effect vs Delta, Plug-in",
       x = "Delta (log scale)",
       y = expression(psi[hat]),
       fill = "Interval Type",
       color = "") +
  theme_minimal()
