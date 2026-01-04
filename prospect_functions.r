library(DEoptim)

# Prospect Theory valuation helper
pt_value <- function(x, a_gain, b_loss, lambda) {
  ifelse(x >= 0,
         x^a_gain,
         -lambda * ((-x)^b_loss))
}

# 1) Objective: Negative log-likelihood (Prospect Theory valued x)
objective_exp_discount_nll_PT <- function(par, data_ppn) {

  data_ppn <- data_ppn[order(data_ppn$TrialNumber), ]

  y <- as.matrix(data_ppn[, c("PAscore", "NAscore")])  # T x 2
  x_raw <- as.numeric(data_ppn$Won)                    # T x 1
  Tn <- nrow(y)

  # unpack core parameters
  alpha <- par[1:2]                      # length 2
  gamma_pa <- par[3]
  gamma_na <- par[4]
  Gamma <- diag(c(gamma_pa, gamma_na))   # 2x2 diagonal
  beta <- matrix(par[5:6], 2, 1)         # 2x1

  # Prospect Theory parameters
  a_gain <- par[7]
  b_loss <- par[8]
  lambda <- par[9]

  # apply PT valuation to x BEFORE recursion
  x <- pt_value(x_raw, a_gain, b_loss, lambda)

  # construct Sigma from g (Cholesky factor)
  g <- matrix(0, 2, 2)
  g[lower.tri(g, diag = TRUE)] <- par[10:12]
  Sigma <- g %*% t(g)

  # discounted sum
  h <- c(0, 0)  # state for discounted accumulation (2-vector PA/NA)

  # precomputing constants for NLL
  d <- 2
  detS <- det(Sigma)
  invS <- solve(Sigma)
  const <- d * log(2*pi) + log(detS)


  nll <- 0
  for (t in 1:Tn) {
    h  <- as.numeric(beta * x[t] + Gamma %*% h)
    mu <- alpha + h
    e  <- y[t, ] - mu
    nll <- nll + 0.5 * (const + as.numeric(t(e) %*% invS %*% e))
  }

  as.numeric(nll)
}

# 2) Objective: Least Squares (Prospect Theory valued x)
objective_exp_discount_LS_PT <- function(par, data_ppn) {

  data_ppn <- data_ppn[order(data_ppn$TrialNumber), ]

  y <- as.matrix(data_ppn[, c("PAscore", "NAscore")])  # T x 2
  x_raw <- as.numeric(data_ppn$Won)
  Tn <- nrow(y)

  # unpack core parameters
  alpha <- par[1:2]
  gamma_pa <- par[3]
  gamma_na <- par[4]
  Gamma <- diag(c(gamma_pa, gamma_na))
  beta <- matrix(par[5:6], 2, 1)

  # Prospect Theory parameters
  a_gain <- par[7]
  b_loss <- par[8]
  lambda <- par[9]

  # apply PT valuation
  x <- pt_value(x_raw, a_gain, b_loss, lambda)

  # recursion
  h <- c(0, 0)
  ssr <- 0

  for (t in 1:Tn) {
    h  <- as.numeric(beta * x[t] + Gamma %*% h)
    mu <- alpha + h
    e  <- y[t, ] - mu
    ssr <- ssr + sum(e^2)
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

# Simulation under Exp Discounting + PT-valued x

simulate_exp_discount_PT <- function(Tn, x_raw, par_true) {

  alpha <- par_true[1:2]
  Gamma <- diag(par_true[3:4])
  beta  <- matrix(par_true[5:6], 2, 1)

  a_gain <- par_true[7]
  b_loss <- par_true[8]
  lambda <- par_true[9]

  # apply PT valuation BEFORE recursion
  x <- pt_value(x_raw, a_gain, b_loss, lambda)

  # noise covariance
  g <- matrix(0, 2, 2)
  g[lower.tri(g, diag = TRUE)] <- par_true[10:12]
  Sigma <- g %*% t(g)

  eps <- mvrnorm(n = Tn, mu = c(0, 0), Sigma = Sigma)

  y <- matrix(NA, nrow = Tn, ncol = 2)
  h <- c(0, 0)

  for (t in 1:Tn) {
    h  <- as.numeric(beta * x[t] + Gamma %*% h)
    mu <- alpha + h
    y[t, ] <- mu + eps[t, ]
  }

  out <- data.frame(
    TrialNumber = 1:Tn,
    Won = x_raw,          # keep raw x in the dataset (objective applies PT)
    PAscore = y[, 1],
    NAscore = y[, 2]
  )
  out
}

# Bounds + parameter sampler
# Parameter order:
# 1:2 alpha (PA, NA)
# 3:4 gamma (PA, NA)
# 5:6 beta  (PA, NA)
# 7:8 r_gain, r_loss
# 9 lambda
# 10:12 g11, g21, g22 (lower-triangular entries of g)

lower <- c(
  -5, -5,          # alpha
   0,  0,          # gamma
  -5, -5,          # beta
   0.10, 0.10,     # r_gain, r_loss
   0.10,           # lambda
   0.05, -2, 0.05  # g11>0, g21 free, g22>0
)

upper <- c(
   5,  5,
  0.999, 0.999,
   5,  5,
   2.00, 2.00,     # r_gain, r_loss
   5.00,           # lambda
   3.00, 2.00, 3.00
)

sample_par <- function(lower, upper) runif(length(lower), min = lower, max = upper)


# DEoptim Estimation with nll

estimate_par_DEoptim_PT <- function(data_sim, lower, upper, itermax = 1500) {
  fit <- DEoptim(
    fn = objective_exp_discount_nll_PT,
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

  # same T and same raw x-series as the template participant
  data_sim <- simulate_exp_discount_PT(Tn = Tn_template, x_raw = x_template, par_true = par_true)

  par_est <- estimate_par_DEoptim_PT(data_sim, lower, upper, itermax = 300)

  true_mat[i, ] <- par_true
  est_mat[i, ]  <- par_est

  cat("Done", i, "of", N, "\n")
}

colnames(true_mat) <- colnames(est_mat) <- c(
  "alpha_PA","alpha_NA",
  "gamma_PA","gamma_NA",
  "beta_PA","beta_NA",
  "a_gain","b_loss","lambda",
  "g11","g21","g22"
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
    title = "Parameter recovery (PT valuation): True vs Estimated (N=100)",
    x = "True value",
    y = "Estimated value"
  ) +
  theme_minimal()
