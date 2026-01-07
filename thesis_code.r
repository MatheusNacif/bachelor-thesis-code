# Libraries
library(DEoptim)
 
# Dataset
data = read.csv('Data_CIAC.csv')

##############################################
# Objective Functions for DEoptim
##############################################
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

  y <- as.matrix(data_ppn[, c("PAscore", "NAscore")]) # T x 2
  x <- as.numeric(data_ppn$Won) # T x 1
  Tn <- nrow(y)

  # necessaary parameters
  alpha <- par[1:2] # length 2
  gamma_pa <- par[3]
  gamma_na <- par[4]
  Gamma <- diag(c(gamma_pa, gamma_na)) # 2x2 diagonal
  beta <- matrix(par[5:6], 2, 1) # 2x1

  # construct Sigma from g (Cholesky factor)
  g <- matrix(0, 2, 2)
  g[lower.tri(g, diag = TRUE)] <- par[7:9]
  Sigma <- g %*% t(g)

  # discounted sum
  h <- c(0, 0) # state for discounted accumulation (2-vector PA/NA)

  # precomputing constants for NLL
  d <- 2
  detS <- det(Sigma)
  invS <- solve(Sigma)
  const <- d * log(2 * pi) + log(detS)

  nll <- 0
  for (t in 1:Tn) {
    # one trial per participant
    h <- as.numeric(beta * x[t] + Gamma %*% h) #updates discounted sum
    mu <- alpha + h # mean vector = baseline + discounted sum
    e <- y[t, ] - mu # resilduals = observed - predicted
    nll <- nll + 0.5 * (const + as.numeric(t(e) %*% invS %*% e)) # accumulate NLL using Mahalanobis distance
  }

  as.numeric(nll)
}

# Addition of Prospect Theory valuation

# Prospect Theory valuation helper
pt_value <- function(x, a_gain, b_loss, lambda) {
  ifelse(x >= 0, abs(x)^a_gain, -lambda * (abs(x)^b_loss))
}

# 2) Objective: Negative log-likelihood (Prospect Theory valued x)
objective_exp_discount_nll_PT <- function(par, data_ppn) {
  data_ppn <- data_ppn[order(data_ppn$TrialNumber), ]

  y <- as.matrix(data_ppn[, c("PAscore", "NAscore")]) # T x 2
  x_raw <- as.numeric(data_ppn$Won) # T x 1
  Tn <- nrow(y)

  # unpack core parameters
  alpha <- par[1:2] # length 2
  gamma_pa <- par[3]
  gamma_na <- par[4]
  Gamma <- diag(c(gamma_pa, gamma_na)) # 2x2 diagonal
  beta <- matrix(par[5:6], 2, 1) # 2x1

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
  h <- c(0, 0) # state for discounted accumulation (2-vector PA/NA)

  # precomputing constants for NLL
  d <- 2
  detS <- det(Sigma)
  invS <- solve(Sigma)
  const <- d * log(2 * pi) + log(detS)

  nll <- 0
  for (t in 1:Tn) {
    h <- as.numeric(beta * x[t] + Gamma %*% h)
    mu <- alpha + h
    e <- y[t, ] - mu
    nll <- nll + 0.5 * (const + as.numeric(t(e) %*% invS %*% e))
  }

  as.numeric(nll)
}

##############################################
# Helper Functions: Fit Metrics and residuals
##############################################

fit_metrics <- function(best_nll, Tn, k) {
  LL  <- -best_nll
  AIC <- 2*k - 2*LL
  BIC <- log(Tn)*k - 2*LL
  c(LL = LL, AIC = AIC, BIC = BIC)
}

# Compute residuals e_t = y_t - mu_t given parameters
residuals_linear <- function(par, data_ppn) {
  data_ppn <- data_ppn[order(data_ppn$TrialNumber), ]
  y <- as.matrix(data_ppn[, c("PAscore", "NAscore")])
  x <- as.numeric(data_ppn$Won)
  Tn <- nrow(y)

  alpha <- par[1:2]
  Gamma <- diag(c(par[3], par[4]))
  beta  <- matrix(par[5:6], 2, 1)

  h <- c(0, 0)
  mu <- matrix(NA_real_, Tn, 2)
  for (t in 1:Tn) {
    h <- as.numeric(beta * x[t] + Gamma %*% h)
    mu[t, ] <- alpha + h
  }
  y - mu
}

residuals_PT <- function(par, data_ppn) {
  data_ppn <- data_ppn[order(data_ppn$TrialNumber), ]
  y <- as.matrix(data_ppn[, c("PAscore", "NAscore")])
  x_raw <- as.numeric(data_ppn$Won)
  Tn <- nrow(y)

  alpha <- par[1:2]
  Gamma <- diag(c(par[3], par[4]))
  beta  <- matrix(par[5:6], 2, 1)

  a_gain <- par[7]
  b_loss <- par[8]
  lambda <- par[9]
  x <- pt_value(x_raw, a_gain, b_loss, lambda)

  h <- c(0, 0)
  mu <- matrix(NA_real_, Tn, 2)
  for (t in 1:Tn) {
    h <- as.numeric(beta * x[t] + Gamma %*% h)
    mu[t, ] <- alpha + h
  }
  y - mu
}

# independence diagnostic: lag-1 autocorrelation for PA and NA_res residuals
lag1_acf <- function(e_mat) {
  PA <- e_mat[,1]; NA_res <- e_mat[,2]
  c(
    acf_PA_lag1 = cor(PA[-1], PA[-length(PA)], use="pairwise.complete.obs"),
    acf_NA_lag1 = cor(NA_res[-1], NA_res[-length(NA_res)], use="pairwise.complete.obs"),
    corr_PA_NA  = cor(PA, NA_res, use="pairwise.complete.obs")
  )
}


##############################################
# DEoptim Bounds
##############################################

# The following bounds are based on prior literature and practical considerations

bounds_linear <- list(
  lower = c(
    -5, -5,        # alpha
     0,  0,        # gamma
    -5, -5,        # beta
     0.05, -2, 0.05  # g11>0, g21 free, g22>0
  ),
  upper = c(
     5,  5,        # alpha
    0.999, 0.999,  # gamma
     5,  5,        # beta
     3,   2,  3    # g entries
  )
)

bounds_pt <- list(
  lower = c(
    -5, -5,          # alpha
     0,  0,          # gamma
    -5, -5,          # beta
     0.001, 0.001,     # r_gain, r_loss
     0.001,           # lambda
     0.05, -2, 0.05  # g11>0, g21 free, g22>0
  ),
  upper = c(
     5,  5,
    0.999, 0.999,
     5,  5,
     1, 1,     # r_gain, r_loss
     3.00,           # lambda
     3,   2,  3      # g entries
  )
)

k_linear <- 9
k_pt     <- 12

# the following is used to set DEoptim population size
# according to their recommendation: NP = 10 * D
p_vec_linear <- length(bounds_linear$lower)
p_vec_pt     <- length(bounds_pt$lower)

##############################################
# Fitting the models to one participant
##############################################

DE_ctrl_stationary_lin <- DEoptim.control(NP = p_vec_linear * 10, itermax = 200, trace = FALSE)
DE_ctrl_stationary_pt <- DEoptim.control(NP = p_vec_pt * 10, itermax = 200, trace = FALSE)

# The following function fits both models to one participant's data
# and returns parameter estimates, fit metrics, LR test, and residual diagnostics

fit_stationary <- function(data_ppn) {
  data_ppn <- data_ppn[order(data_ppn$TrialNumber), ]
  Tn <- nrow(data_ppn)

  fit_lin <- DEoptim(
    fn = objective_exp_discount_nll,
    lower = bounds_linear$lower,
    upper = bounds_linear$upper,
    data_ppn = data_ppn,
    control = DE_ctrl_stationary_lin
  )
  par_lin <- fit_lin$optim$bestmem
  nll_lin <- fit_lin$optim$bestval
  met_lin <- fit_metrics(nll_lin, Tn, k_linear)

  fit_ptm <- DEoptim(
    fn = objective_exp_discount_nll_PT,
    lower = bounds_pt$lower,
    upper = bounds_pt$upper,
    data_ppn = data_ppn,
    control = DE_ctrl_stationary_pt
  )
  par_pt <- fit_ptm$optim$bestmem
  nll_pt <- fit_ptm$optim$bestval
  met_pt <- fit_metrics(nll_pt, Tn, k_pt)

  # Likelihood ratio test (nested), df = 3
  LR   <- 2 * (met_pt["LL"] - met_lin["LL"])
  p_LR <- pchisq(LR, df = 3, lower.tail = FALSE)

  # Residual diagnostics
  diag_lin <- lag1_acf(residuals_linear(par_lin, data_ppn))
  diag_pt  <- lag1_acf(residuals_PT(par_pt, data_ppn))

  list(
    par_lin = par_lin, metrics_lin = met_lin,
    par_pt  = par_pt,  metrics_pt  = met_pt,
    LR = as.numeric(LR), p_LR = as.numeric(p_LR),
    resid_diag_lin = diag_lin,
    resid_diag_pt  = diag_pt
  )
}

###############################################
# Moving Window Estimation
###############################################

DE_ctrl_window_lin <- DEoptim.control(NP = p_vec_linear * 10, itermax = 120, trace = FALSE)
DE_ctrl_window_pt <- DEoptim.control(NP = p_vec_pt * 10, itermax = 120, trace = FALSE)

# The following function performs moving-window estimation for one participant
# and returns a data frame with parameter trajectories and nll per window

moving_window <- function(data_ppn, model = c("linear","pt"), window = 60) {
  model <- match.arg(model)
  data_ppn <- data_ppn[order(data_ppn$TrialNumber), ]
  Tn <- nrow(data_ppn)
  if (Tn < window) return(NULL)

  nW <- Tn - window + 1

  if (model == "linear") {
    traj <- matrix(NA_real_, nW, k_linear)
    nll  <- rep(NA_real_, nW)
    for (w in 1:nW) {
      seg <- data_ppn[w:(w+window-1), ]
      fit <- DEoptim(objective_exp_discount_nll,
                     bounds_linear$lower, bounds_linear$upper,
                     data_ppn = seg, control = DE_ctrl_window_lin)
      traj[w,] <- fit$optim$bestmem
      nll[w]   <- fit$optim$bestval
    }
    colnames(traj) <- c("alpha_PA","alpha_NA","gamma_PA","gamma_NA","beta_PA","beta_NA","g11","g21","g22")
  } else {
    traj <- matrix(NA_real_, nW, k_pt)
    nll  <- rep(NA_real_, nW)
    for (w in 1:nW) {
      seg <- data_ppn[w:(w+window-1), ]
      fit <- DEoptim(objective_exp_discount_nll_PT,
                     bounds_pt$lower, bounds_pt$upper,
                     data_ppn = seg, control = DE_ctrl_window_pt)
      traj[w,] <- fit$optim$bestmem
      nll[w]   <- fit$optim$bestval
    }
    colnames(traj) <- c("alpha_PA","alpha_NA","gamma_PA","gamma_NA","beta_PA","beta_NA",
                        "r_plus","r_minus","lambda","g11","g21","g22")
  }

  data.frame(
    start_trial = data_ppn$TrialNumber[1:nW],
    end_trial   = data_ppn$TrialNumber[window:(window+nW-1)],
    traj,
    nll = nll
  )
}

################################################
# Full Analysis Across Participants
################################################

# The following function runs the full analysis across all participants
# It fits both models, computes fit metrics, LR tests, residual diagnostics,
# and moving-window trajectories for each participant

run_analysis <- function(data, id_col = "Ppn", window = 70) {
  ids <- unique(data[[id_col]])

  summary <- data.frame(
    Ppn = ids,
    LL_lin = NA_real_, AIC_lin = NA_real_, BIC_lin = NA_real_,
    LL_pt  = NA_real_, AIC_pt  = NA_real_, BIC_pt  = NA_real_,
    LR = NA_real_, p_LR = NA_real_,
    acfPA_lin = NA_real_, acfNA_lin = NA_real_, corr_lin = NA_real_,
    acfPA_pt  = NA_real_, acfNA_pt  = NA_real_, corr_pt  = NA_real_
  )

  window_lin <- vector("list", length(ids)); names(window_lin) <- as.character(ids)
  window_pt  <- vector("list", length(ids)); names(window_pt)  <- as.character(ids)

  fits <- vector("list", length(ids)); names(fits) <- as.character(ids)

  for (i in seq_along(ids)) {
    pid <- ids[i]
    data_ppn <- subset(data, data[[id_col]] == pid)

    cat("Fitting Ppn:", pid, "\n")
    out <- fit_stationary(data_ppn)
    fits[[as.character(pid)]] <- out

    summary$LL_lin[i]  <- out$metrics_lin["LL"]
    summary$AIC_lin[i] <- out$metrics_lin["AIC"]
    summary$BIC_lin[i] <- out$metrics_lin["BIC"]

    summary$LL_pt[i]   <- out$metrics_pt["LL"]
    summary$AIC_pt[i]  <- out$metrics_pt["AIC"]
    summary$BIC_pt[i]  <- out$metrics_pt["BIC"]

    summary$LR[i]   <- out$LR
    summary$p_LR[i] <- out$p_LR

    summary$acfPA_lin[i] <- out$resid_diag_lin["acf_PA_lag1"]
    summary$acfNA_lin[i] <- out$resid_diag_lin["acf_NA_lag1"]
    summary$corr_lin[i]  <- out$resid_diag_lin["corr_PA_NA"]

    summary$acfPA_pt[i]  <- out$resid_diag_pt["acf_PA_lag1"]
    summary$acfNA_pt[i]  <- out$resid_diag_pt["acf_NA_lag1"]
    summary$corr_pt[i]   <- out$resid_diag_pt["corr_PA_NA"]

    # Moving window trajectories
    window_lin[[as.character(pid)]] <- moving_window(data_ppn, "linear", window)
    window_pt[[as.character(pid)]]  <- moving_window(data_ppn, "pt", window)
  }

  summary$PT_better_AIC <- summary$AIC_pt < summary$AIC_lin
  summary$PT_better_BIC <- summary$BIC_pt < summary$BIC_lin

  list(summary = summary, fits = fits, mw_linear = window_lin, mw_pt = window_pt)
}

##############################
# Running it
###############################

data <- read.csv("Data_CIAC.csv")
results <- run_analysis(data, id_col = "Ppn", window = 70)

results$summary
write.csv(results$summary, "model_comparison_summary.csv", row.names = FALSE)

##############################
# Analysis
###############################

# Does PT improve explanatory power under stationarity?

mean(results$summary$AIC_pt < results$summary$AIC_lin)
mean(results$summary$BIC_pt < results$summary$BIC_lin)

# Percentage of participants with significant LR test for PT over linear

mean(results$summary$p_LR < .05)

# Residual Diagnostics for model fit

summary(results$summary$acfPA_lin - results$summary$acfPA_pt)
summary(results$summary$acfNA_lin - results$summary$acfNA_pt)

# Temporal Instability Analysis
# # Plots 3 participants: highest variability, lowest variability, and average (median) variability
params_to_check <- c("beta_PA", "beta_NA", "gamma_PA", "gamma_NA")

# Reference parameter for ranking variability
ref_param <- "beta_PA"

# Compute within-participant SD for the reference parameter (PT model used for ranking)
sd_ref <- sapply(results$mw_pt, function(x) {
  if (!is.null(x) && ref_param %in% names(x)) sd(x[[ref_param]], na.rm = TRUE) else NA_real_
})
sd_ref <- sd_ref[!is.na(sd_ref)]
sd_ref_sorted <- sort(sd_ref)

# Select participants
pid_low  <- names(sd_ref_sorted)[1]
pid_high <- names(sd_ref_sorted)[length(sd_ref_sorted)]
pid_mid  <- names(sd_ref_sorted)[ceiling(length(sd_ref_sorted) / 2)]

# Fixed order for plotting and legend
selected_pids <- c(pid_low, pid_mid, pid_high)

# Color mapping by variability group
group_cols <- c(
  "Lowest variability"  = "black",
  "Average variability" = "green3",
  "Highest variability" = "red"
)

# Plot settings: fixed x-axis from known window count
window_size <- 60
Tn_total <- 152
n_windows <- Tn_total - window_size + 1  # 93
x_rng <- c(1, n_windows)

# Helper: adds selected participants with fixed colors
add_selected_lines_by_group <- function(mw_list, param, selected_pids, group_cols, lwd = 2) {
  for (i in seq_along(selected_pids)) {
    pid <- selected_pids[i]
    mw <- mw_list[[pid]]
    if (!is.null(mw) && param %in% names(mw)) {
      col <- unname(group_cols[i])
      lines(mw$start_trial, mw[[param]], col = col, lwd = lwd)
    }
  }
}

# Plotting: two plots per parameter (Linear then PT), pooled y-range for comparability
for (param in params_to_check) {

  # pooled y-range across models, for selected participants
  vals_lin <- unlist(lapply(selected_pids, function(pid) {
    mw <- results$mw_linear[[pid]]
    if (!is.null(mw)) mw[[param]] else NA_real_
  }))
  vals_pt <- unlist(lapply(selected_pids, function(pid) {
    mw <- results$mw_pt[[pid]]
    if (!is.null(mw)) mw[[param]] else NA_real_
  }))
  y_rng <- range(c(vals_lin, vals_pt), na.rm = TRUE)

  # ---- Linear plot
  plot(
    NA,
    xlim = x_rng,
    ylim = y_rng,
    xlab = "Window start trial",
    ylab = param,
    main = paste("Moving-window trajectories (Linear):", param)
  )

  add_selected_lines_by_group(
    mw_list = results$mw_linear,
    param = param,
    selected_pids = selected_pids,
    group_cols = group_cols
  )

  legend(
    "topright",
    legend = names(group_cols),
    col = group_cols,
    lty = 1,
    lwd = 2,
    cex = 0.9,
    bty = "n"
  )

  # ---- PT plot
  plot(
    NA,
    xlim = x_rng,
    ylim = y_rng,
    xlab = "Window start trial",
    ylab = param,
    main = paste("Moving-window trajectories (Prospect Theory):", param)
  )

  add_selected_lines_by_group(
    mw_list = results$mw_pt,
    param = param,
    selected_pids = selected_pids,
    group_cols = group_cols
  )

  legend(
    "topright",
    legend = names(group_cols),
    col = group_cols,
    lty = 1,
    lwd = 2,
    cex = 0.9,
    bty = "n"
  )
}

# Within-participant variability across windows:
# compute SD per parameter (per participant) for both models

sd_by_param <- data.frame(
  parameter = params_to_check,
  mean_sd_linear = NA_real_,
  mean_sd_pt     = NA_real_,
  prop_pt_more_variable = NA_real_  # proportion of participants where SD_PT > SD_Linear
)

for (k in seq_along(params_to_check)) {
  param <- params_to_check[k]

  sd_lin <- sapply(results$mw_linear, function(x) if (!is.null(x)) sd(x[[param]], na.rm = TRUE) else NA_real_)
  sd_pt  <- sapply(results$mw_pt,     function(x) if (!is.null(x)) sd(x[[param]], na.rm = TRUE) else NA_real_)

  sd_by_param$mean_sd_linear[k] <- mean(sd_lin, na.rm = TRUE)
  sd_by_param$mean_sd_pt[k]     <- mean(sd_pt,  na.rm = TRUE)
  sd_by_param$prop_pt_more_variable[k] <- mean(sd_pt > sd_lin, na.rm = TRUE)
}

sd_by_param

# End of code