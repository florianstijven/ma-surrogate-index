stan_code_prop_line = 
  "
    data {
      int<lower=1> N;              // number of studies
      vector[2] y[N];              // observed effects: [S, T]
      cov_matrix[2] Sigma[N];      // known within-trial covariance matrices
    }
    parameters {
      real mu;                     // mean treatment effect on surrogate
      vector<lower=0>[2] tau;      // between-trial SDs
      real<lower=-1, upper=1> rho; // trial-level correlation
      vector[2] theta[N];          // latent true effects per trial
    }
    model {
      mu ~ normal(log(0.5), 1);           
      tau ~ normal(0, 2);        
      rho ~ uniform(-1, 1);        

      matrix[2,2] Omega;
      Omega[1,1] = 1;
      Omega[1,2] = rho;
      Omega[2,1] = rho;
      Omega[2,2] = 1;

      real mu_beta;
      mu_beta = rho * (tau[2] / tau[1]) * mu;

      vector[2] mu_full;
      mu_full[1] = mu;
      mu_full[2] = mu_beta;

      for (i in 1:N)
        theta[i] ~ multi_normal(mu_full, quad_form_diag(Omega, tau));

      for (i in 1:N)
        y[i] ~ multi_normal(theta[i], Sigma[i]);
    }
    "
stan_code_default = 
  "
    data {
      int<lower=1> N;              // number of studies
      vector[2] y[N];              // observed effects: [S, T]
      cov_matrix[2] Sigma[N];      // known within-trial covariance matrices
    }
    parameters {
      vector[2] mu;                     // mean treatment effect on surrogate and clinical endpoint
      vector<lower=0>[2] tau;      // between-trial SDs
      real<lower=-1, upper=1> rho; // trial-level correlation
      vector[2] theta[N];          // latent true effects per trial
    }
    model {
      mu ~ normal(log(0.5), 1);           
      tau ~ normal(0, 2);        
      rho ~ uniform(-1, 1);        

      matrix[2,2] Omega;
      Omega[1,1] = 1;
      Omega[1,2] = rho;
      Omega[2,1] = rho;
      Omega[2,2] = 1;

      for (i in 1:N)
        theta[i] ~ multi_normal(mu, quad_form_diag(Omega, tau));

      for (i in 1:N)
        y[i] ~ multi_normal(theta[i], Sigma[i]);
    }
    "

stan_model_prop_line = stan_model(model_code = stan_code_prop_line)
stan_model_default = stan_model(model_code = stan_code_default)

fit_surrogacy_model <- function(data, assume_proportional_line, iter = 10000, warmup = 5000, chains = 4, seed = 123) {
  # Define Stan code dynamically based on the assumption of an identity regression line
  stan_model <- if (assume_proportional_line) {
    stan_model_prop_line
  } else {
    stan_model_default
  }
  
  # Prepare input data for Stan
  N <- nrow(data)
  y_matrix <- cbind(data$trt_effect_surrogate_index_est, data$log_RR_est)
  
  Sigma_array <- array(NA, dim = c(N, 2, 2))
  for (i in 1:N) {
    Sigma_array[i,,] <- as.matrix(Matrix::nearPD(data$vcov[[i]])$mat)
  }
  
  stan_data <- list(
    N = N,
    y = y_matrix,
    Sigma = Sigma_array
  )
  
  # Fit model
  fit <- sampling(
    object = stan_model,
    data = stan_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    seed = seed, refresh = FALSE
  )
  
  return(fit)
}
