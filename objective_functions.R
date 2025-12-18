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