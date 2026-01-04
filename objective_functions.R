library(DEoptim)

# This is the example used in the workshop
# y = beta_0 + beta_1 * x + error

# Least squares objective function
objective_function_ls <- function(parameters, data) {
  # Predicted values
  y_pred <- parameters[1] + parameters[2] * data$x
  
  # Residuals
  residuals <- data$y - y_pred
  
  # Sum of squared residuals
  ssr <- sum(residuals^2)
  
  return(ssr)
}

# Maximum likelihood objective function assuming normal errors
objective_function_mll <- function(parameters, data) {
  # Predicted values
  y_pred <- parameters[1] + parameters[2] * data$x
  
  # Get likelihoods for each datapoint
  L <- dnorm(
    data$y, 
    mean = y_pred,
    sd = parameters[3]
  )
  
  # Transform to min-log-likelihood
  return(-sum(log(L)))
}

data = read.csv('Data_CIAC.csv')

#    Model: Exponential Discounting model (multi-dimensional output)
#      y_t = alpha + sum_{j>=0} Gamma^j beta x_{t-j} + eps_t
#    with recursion:
#      h_t = beta*x_t + Gamma*h_{t-1}
#      mu_t = alpha + h_t
#
#    d = 2 (PA, NA)
#    m = 1 (Won)
#    Gamma is restricted to diagonal (two forgetting factors)

# 1) Objective: Negative log-likelihood

objective_exp_discount_nll <- function(par, data_ppn) {

  data_ppn <- data_ppn[order(data_ppn$TrialNumber), ]

  y <- as.matrix(data_ppn[, c("PAscore", "NAscore")])  # T x 2
  x <- as.numeric(data_ppn$Won)                        # T x 1
  Tn <- nrow(y)

  # necessaary parameters
  alpha <- par[1:2]                    # length 2
  gamma_pa <- par[3]
  gamma_na <- par[4]
  Gamma <- diag(c(gamma_pa, gamma_na)) # 2x2 diagonal
  beta <- matrix(par[5:6], 2, 1)  # 2x1

  # construct Sigma from g (Cholesky factor)
  g <- matrix(0, 2, 2)
  g[lower.tri(g, diag = TRUE)] <- par[7:9]
  Sigma <- g %*% t(g)

  # discounted sum
  h <- c(0, 0)  # state for discounted accumulation (2-vector PA/NA)

  # precomputing constants for NLL
  d <- 2
  detS <- det(Sigma)
  invS <- solve(Sigma)
  const <- d * log(2*pi) + log(detS)

  nll <- 0
  for (t in 1:Tn) {        # one trial per participant
    h <- as.numeric(beta * x[t] + Gamma %*% h)    #updates discounted sum
    mu <- alpha + h                          # mean vector = baseline + discounted sum
    e  <- y[t, ] - mu                        # resilduals = observed - predicted
    nll <- nll + 0.5 * (const + as.numeric(t(e) %*% invS %*% e))  # accumulate NLL using Mahalanobis distance
  }

  as.numeric(nll)
}

# 2) Objective: Least Squares

objective_exp_discount_LS <- function(par, data_ppn) {

  # Order trials
  data_ppn <- data_ppn[order(data_ppn$TrialNumber), ]

  y <- as.matrix(data_ppn[, c("PAscore", "NAscore")])  # T x 2
  x <- as.numeric(data_ppn$Won)
  Tn <- nrow(y)

  # unpack parameters
  alpha <- par[1:2]                    # length 2
  gamma_pa <- par[3]
  gamma_na <- par[4]
  Gamma <- diag(c(gamma_pa, gamma_na)) # 2x2 diagonal
  beta <- matrix(par[5:6], 2,  1)  # 2x1

  # recursion for discounted sum
  h <- c(0, 0)
  ssr <- 0

  for (t in 1:Tn) {            # one trial per participant
    h <- as.numeric(beta * x[t] + Gamma %*% h)   # update discounted sum
    mu <- alpha + h                         # mean vector = baseline + discounted sum
    e  <- y[t, ] - mu                       # resilduals = observed - predicted
    ssr <- ssr + sum(e^2)                   # accumulate sum of squared residuals to be minimized
  }

  as.numeric(ssr)
}

# Recovery study with one participant
library(MASS)
library(ggplot2)

data = read.csv('Data_CIAC.csv')

ppn_id <- data$Ppn[1]
data_ppn_template <- subset(data, Ppn == ppn_id)

Tn_template <- nrow(data_ppn_template)
x_template  <- as.numeric(data_ppn_template$Won)

# Simulation under exponential discounting model

simulate_exp_discount <- function(Tn, x, par_true) {

  alpha <- par_true[1:2]
  Gamma <- diag(par_true[3:4])
  beta  <- matrix(par_true[5:6], 2, 1)

  g <- matrix(0, 2, 2)
  g[lower.tri(g, diag = TRUE)] <- par_true[7:9]
  Sigma <- g %*% t(g)

  # simulate eps_t ~ MVN(0, Sigma)
  eps <- mvrnorm(n = Tn, mu = c(0,0), Sigma = Sigma)

  y <- matrix(NA, nrow = Tn, ncol = 2)
  h <- c(0, 0)

  for (t in 1:Tn) {
    h <- as.numeric(beta * x[t] + Gamma %*% h)
    mu <- alpha + h
    y[t, ] <- mu + eps[t, ]
  }

  # return a data frame like the original data
  out <- data.frame(
    TrialNumber = 1:Tn,
    Won = x,
    PAscore = y[,1],
    NAscore = y[,2]
  )
  out
}

# Bounds + parameter sampler
# Parameter order:
# 1:2 alpha (PA, NA)
# 3:4 gamma (PA, NA)
# 5:6 beta  (PA, NA)
# 7:9 g11, g21, g22 (lower-triangular entries of g)

lower <- c(
  -5, -5,        # alpha
   0,  0,        # gamma
  -5, -5,        # beta
   0.05, -2, 0.05  # g11>0, g21 free, g22>0
)

upper <- c(
   5,  5,        # alpha
  0.999, 0.999,  # gamma
   5,  5,        # beta
   3,   2,  3    # g entries
)

sample_par <- function(lower, upper) runif(length(lower), min = lower, max = upper)


# DEoptim Estimation with nll

estimate_par_DEoptim <- function(data_sim, lower, upper, itermax = 1500) {
  fit <- DEoptim(
    fn = objective_exp_discount_nll,
    lower = lower,
    upper = upper,
    control = DEoptim.control(itermax = itermax, trace = FALSE),
    data_ppn = data_sim
  )
  fit$optim$bestmem
}


# Recovery study loop (N = 100)

set.seed(123)

N <- 100
true_mat <- matrix(NA, nrow = N, ncol = length(lower))
est_mat  <- matrix(NA, nrow = N, ncol = length(lower))

for (i in 1:N) {
  par_true <- sample_par(lower, upper)

  # Use the same T and the same x series as the study participant
  data_sim <- simulate_exp_discount(Tn = Tn_template, x = x_template, par_true = par_true)

  par_est <- estimate_par_DEoptim(data_sim, lower, upper, itermax = 300)

  true_mat[i, ] <- par_true
  est_mat[i, ]  <- par_est

  cat("Done", i, "of", N, "\n")
}

colnames(true_mat) <- colnames(est_mat) <- c(
  "alpha_PA","alpha_NA","gamma_PA","gamma_NA","beta_PA","beta_NA","g11","g21","g22"
)

df_rec <- data.frame(
  rep = 1:N,
  true_mat,
  est_mat
)

# Long format for plotting
df_long <- do.call(rbind, lapply(colnames(true_mat), function(p) {
  data.frame(
    param = p,
    true  = true_mat[, p],
    est   = est_mat[, p]
  )
}))

# Visualization: true vs estimated with y=x line

ggplot(df_long, aes(x = true, y = est)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ param, scales = "free") +
  labs(
    title = "Parameter recovery: True vs Estimated (N=100)",
    x = "True value",
    y = "Estimated value"
  ) +
  theme_minimal()
