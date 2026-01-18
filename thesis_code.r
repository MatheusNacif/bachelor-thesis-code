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

  # Cholesky factorization of Sigma (fail-fast if invalid)
  R <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(R)) return(1e12)  # heavy penalty keeps DEoptim away

  # log(det(Sigma)) via Cholesky
  logdetS <- 2 * sum(log(diag(R)))

  # constants (d = 2)
  const <- 2 * log(2 * pi) + logdetS

  # discounted sum state
  h <- c(0, 0)

  nll <- 0
  for (t in 1:Tn) {
    h  <- as.numeric(beta * x[t] + Gamma %*% h)
    mu <- alpha + h
    e  <- y[t, ] - mu

    # quadratic form e' Sigma^{-1} e via Cholesky:
    # if Sigma = R'R, then solve(Sigma, e) done through backsolve
    z <- backsolve(R, e, transpose = TRUE)
    qf <- sum(z * z)

    nll <- nll + 0.5 * (const + qf)
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

  # Cholesky factorization of Sigma (fail-fast if invalid)
  R <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(R)) return(1e12)

  # log(det(Sigma)) via Cholesky
  logdetS <- 2 * sum(log(diag(R)))

  # constants (d = 2)
  const <- 2 * log(2 * pi) + logdetS

  # discounted sum state
  h <- c(0, 0)

  nll <- 0
  for (t in 1:Tn) {
    h  <- as.numeric(beta * x[t] + Gamma %*% h)
    mu <- alpha + h
    e  <- y[t, ] - mu

    # quadratic form via Cholesky
    z <- backsolve(R, e, transpose = TRUE)
    qf <- sum(z * z)

    nll <- nll + 0.5 * (const + qf)
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
    0, 0,        # alpha
     0,  0,        # gamma
    -0.333, -0.333,        # beta
     1e-5, -1, 1e-5  # g11>0, g21 free, g22>0
  ),
  upper = c(
    1,  1,        # alpha
      0.999, 0.999,  # gamma
    0.333,  0.333,        # beta
     0.5,   1,  0.5    # g entries
  )
)

bounds_pt <- list(
  lower = c(
    0, 0,          # alpha
     0,  0,          # gamma
    -0.333, -0.333,          # beta
     0.001, 0.001,     # r_gain, r_loss
    0.001,           # lambda
     1e-5, -1, 1e-5  # g11>0, g21 free, g22>0
  ),
  upper = c(
    1,  1,
      0.999, 0.999,
    0.333,  0.333,
     1, 1,     # r_gain, r_loss
    3.00,           # lambda
     0.5,   1,  0.5      # g entries
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

DE_ctrl_stationary_lin <- DEoptim.control(strategy = 6, NP = p_vec_pt * 20, itermax = 500, trace = FALSE, p = 0.5, steptol = 50)
DE_ctrl_stationary_pt <- DEoptim.control(strategy = 6, NP = p_vec_pt * 20, itermax = 500, trace = FALSE, p = 0.5, steptol = 50)

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

DE_ctrl_window_lin <- DEoptim.control(strategy = 6, NP = p_vec_pt * 10, itermax = 150, trace = FALSE, p = 0.5, steptol = 25)
DE_ctrl_window_pt <- DEoptim.control(strategy = 6, NP = p_vec_pt * 10, itermax = 150, trace = FALSE, p = 0.5, steptol = 25)


# The following function performs moving-window estimation for one participant
# and returns a data frame with parameter trajectories and nll per window.
# The function uses the last best parameters as initialization for the next window.

moving_window <- function(data_ppn, model = c("linear","pt"), window = 60, step = 1) {
  model <- match.arg(model)
  data_ppn <- data_ppn[order(data_ppn$TrialNumber), ]
  Tn <- nrow(data_ppn)
  if (Tn < window) return(NULL)

  # window start indices with step
  starts <- seq(1, Tn - window + 1, by = step)
  nW <- length(starts)

  if (model == "linear") {
    lower <- bounds_linear$lower; upper <- bounds_linear$upper
    ctrl  <- DE_ctrl_window_lin
    k     <- k_linear

    traj <- matrix(NA_real_, nW, k)
    nll  <- rep(NA_real_, nW)

    # initial population (created ONCE)
    D  <- length(lower)
    NP <- ctrl$NP
    initialpop <- matrix(runif(NP * D, min = lower, max = upper), nrow = NP, byrow = TRUE)

    # explicitly randomize row 1
    initialpop[1, ] <- runif(D, min = lower, max = upper)

    for (i in 1:nW) {
      w <- starts[i]
      seg <- data_ppn[w:(w + window - 1), ]

      ctrl$initialpop <- initialpop

      fit <- DEoptim(objective_exp_discount_nll,
                     lower = lower, upper = upper,
                     data_ppn = seg, control = ctrl)

      best <- fit$optim$bestmem

      traj[i, ] <- best
      nll[i]    <- fit$optim$bestval

      # update row 1 for the NEXT window (clip to bounds)
      initialpop[1, ] <- pmin(pmax(best, lower), upper)
    }

    colnames(traj) <- c("alpha_PA","alpha_NA","gamma_PA","gamma_NA","beta_PA","beta_NA","g11","g21","g22")

  } else {
    lower <- bounds_pt$lower; upper <- bounds_pt$upper
    ctrl  <- DE_ctrl_window_pt
    k     <- k_pt

    traj <- matrix(NA_real_, nW, k)
    nll  <- rep(NA_real_, nW)

    # initial population (created ONCE)
    D  <- length(lower)
    NP <- ctrl$NP
    initialpop <- matrix(runif(NP * D, min = lower, max = upper), nrow = NP, byrow = TRUE)

    # explicitly randomize row 1
    initialpop[1, ] <- runif(D, min = lower, max = upper)

    for (i in 1:nW) {
      w <- starts[i]
      seg <- data_ppn[w:(w + window - 1), ]

      ctrl$initialpop <- initialpop

      fit <- DEoptim(objective_exp_discount_nll_PT,
                     lower = lower, upper = upper,
                     data_ppn = seg, control = ctrl)

      best <- fit$optim$bestmem

      traj[i, ] <- best
      nll[i]    <- fit$optim$bestval

      # update row 1 for next window (clip)
      initialpop[1, ] <- pmin(pmax(best, lower), upper)
    }

    colnames(traj) <- c("alpha_PA","alpha_NA","gamma_PA","gamma_NA","beta_PA","beta_NA",
                        "r_plus","r_minus","lambda","g11","g21","g22")
  }

  data.frame(
    start_trial = data_ppn$TrialNumber[starts],
    end_trial   = data_ppn$TrialNumber[starts + window - 1],
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

run_analysis <- function(data, id_col = "Ppn", window = 60) {
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
# Time testing
##############################

#pid_test <- unique(data$Ppn)[1]
#data_test <- subset(data, Ppn == pid_test)

# stationary fits + both moving-window fits
#t0 <- Sys.time()

#out_stationary <- fit_stationary(data_test)
#out_mw_lin     <- moving_window(data_test, model = "linear", window = 70)
#out_mw_pt      <- moving_window(data_test, model = "pt",     window = 70)

#t1 <- Sys.time()

#elapsed_sec <- as.numeric(difftime(t1, t0, units = "secs"))
#cat("Timing test for Ppn =", pid_test, "\n")
#cat("Elapsed time:", round(elapsed_sec, 2), "seconds (", round(elapsed_sec/60, 2), "minutes )\n")

##############################
# Running it
###############################

data <- read.csv("Data_CIAC.csv")
results <- run_analysis(data, id_col = "Ppn", window = 60)

results$summary
write.csv(results$summary, "model_comparison_summary.csv", row.names = FALSE)

##############################
# Descriptive statistics
###############################

# Affect descriptives (trial-level)

affect_desc <- data.frame(
  Measure = c("PA", "NA"),
  Mean = c(mean(data$PAscore, na.rm = TRUE),
           mean(data$NAscore, na.rm = TRUE)),
  SD = c(sd(data$PAscore, na.rm = TRUE),
         sd(data$NAscore, na.rm = TRUE)),
  Min = c(min(data$PAscore, na.rm = TRUE),
          min(data$NAscore, na.rm = TRUE)),
  Max = c(max(data$PAscore, na.rm = TRUE),
          max(data$NAscore, na.rm = TRUE))
)

affect_desc

# Outcome descriptives

outcome_desc <- data.frame(
  Mean = mean(data$Won),
  SD   = sd(data$Won),
  Min  = min(data$Won),
  Max  = max(data$Won),
  Prop_Wins  = mean(data$Won > 0),
  Prop_Losses = mean(data$Won < 0)
)

outcome_desc

# Extract stationary parameter estimates

extract_pars <- function(fits, model = c("linear", "pt")) {
  model <- match.arg(model)

  par_mat <- do.call(rbind, lapply(fits, function(x) {
    if (model == "linear") x$par_lin else x$par_pt
  }))

  as.data.frame(par_mat)
}

pars_linear <- extract_pars(results$fits, "linear")
pars_pt     <- extract_pars(results$fits, "pt")

linear_desc <- data.frame(
  Parameter = colnames(pars_linear),
  Mean = apply(pars_linear, 2, mean, na.rm = TRUE),
  SD   = apply(pars_linear, 2, sd,   na.rm = TRUE),
  Min  = apply(pars_linear, 2, min,  na.rm = TRUE),
  Max  = apply(pars_linear, 2, max,  na.rm = TRUE)
)

linear_desc

pt_desc <- data.frame(
  Parameter = colnames(pars_pt),
  Mean = apply(pars_pt, 2, mean, na.rm = TRUE),
  SD   = apply(pars_pt, 2, sd,   na.rm = TRUE),
  Min  = apply(pars_pt, 2, min,  na.rm = TRUE),
  Max  = apply(pars_pt, 2, max,  na.rm = TRUE)
)

pt_desc

# Descriptives to CSV for Word
write.csv(affect_desc, "descriptives_affect.csv", row.names = FALSE)
write.csv(outcome_desc, "descriptives_outcomes.csv", row.names = FALSE)
write.csv(linear_desc, "descriptives_parameters_linear.csv", row.names = FALSE)
write.csv(pt_desc, "descriptives_parameters_PT.csv", row.names = FALSE)

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

params_to_check <- c("beta_PA", "beta_NA", "gamma_PA", "gamma_NA")

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

# Settings for transfering to Word

theme_word <- theme_bw() +
  theme(
    # Main title & subtitle
    plot.title    = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),

    # Axis titles
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),

    # Axis tick labels
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),

    # Facet strip titles
    strip.text = element_text(size = 14, face = "bold"),

    # Legend
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),

    # Legend spacing
    legend.key.size = unit(1.2, "lines")
  )


# Plot 1: within-participant SD across moving windows

library(ggplot2)

params_to_check <- c("beta_PA", "beta_NA", "gamma_PA", "gamma_NA")

# Compute SD across windows for one parameter, for all participants in a mw_list
get_sd_vector <- function(mw_list, param) {
  sapply(mw_list, function(df) {
    if (is.null(df) || !(param %in% names(df))) return(NA_real_)
    sd(df[[param]], na.rm = TRUE)
  })
}

# One row = one participant x one parameter x one model
sd_long <- do.call(rbind, lapply(params_to_check, function(param) {

  sd_lin <- get_sd_vector(results$mw_linear, param)
  sd_pt  <- get_sd_vector(results$mw_pt,     param)

  # Keep only participants with non-missing SD in both models
  common_ids <- intersect(names(sd_lin)[!is.na(sd_lin)], names(sd_pt)[!is.na(sd_pt)])

  data.frame(
    Ppn = rep(common_ids, times = 2),
    parameter = param,
    model = rep(c("Linear", "Prospect Theory"), each = length(common_ids)),
    sd_within = c(sd_lin[common_ids], sd_pt[common_ids]),
    row.names = NULL
  )
}))

# Make parameter labels prettier for plotting
sd_long$parameter <- factor(sd_long$parameter,
                           levels = params_to_check,
                           labels = c(expression(beta[PA]),
                                      expression(beta["NA"]),
                                      expression(gamma[PA]),
                                      expression(gamma["NA"])))

# Plot
p1 <- ggplot(sd_long, aes(x = model, y = sd_within)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 1.2) +
  facet_wrap(~ parameter, scales = "free_y", labeller = label_parsed) +
  labs(
    x = NULL,
    y = "Within-participant SD across windows",
    title = "Temporal variability of moving-window parameter estimates",
    subtitle = "Each point is one participant (SD computed across all moving windows)"
  ) +
  theme_word +
  theme(legend.position = "none")

print(p1)


# Plot 2:Average moving-window trajectories across participants

mw_to_long_multi <- function(mw_list, model_name, params) {
  do.call(rbind, lapply(names(mw_list), function(pid) {
    df <- mw_list[[pid]]
    if (is.null(df)) return(NULL)

    # keep only params that exist
    keep <- intersect(params, names(df))
    if (length(keep) == 0) return(NULL)

    do.call(rbind, lapply(keep, function(p) {
      data.frame(
        Ppn = pid,
        model = model_name,
        start_trial = df$start_trial,
        parameter = p,
        value = df[[p]],
        row.names = NULL
      )
    }))
  }))
}

df_long <- rbind(
  mw_to_long_multi(results$mw_linear, "Linear", params_to_check),
  mw_to_long_multi(results$mw_pt, "Prospect Theory", params_to_check)
)

# Summarize across participants per (model x start_trial x parameter)
df_sum <- aggregate(
  value ~ model + start_trial + parameter,
  data = df_long,
  FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
)

# Unpack the matrix column
df_sum$mean <- df_sum$value[, "mean"]
df_sum$sd   <- df_sum$value[, "sd"]
df_sum$value <- NULL

# Nicer facet labels
df_sum$parameter <- factor(
  df_sum$parameter,
  levels = params_to_check,
  labels = c(expression(beta[PA]), expression(beta["NA"]),
             expression(gamma[PA]), expression(gamma["NA"]))
)

p2 <- ggplot(df_sum, aes(x = start_trial, y = mean, linetype = model)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.15) +
  geom_line(linewidth = 1) +
  facet_wrap(~ parameter, scales = "free_y", labeller = label_parsed) +
  labs(
    x = "Window start trial",
    y = "Parameter estimate",
    title = "Average moving-window trajectories across participants",
    subtitle = "Line = mean; band = ±1 SD"
  ) +
  theme_word

print(p2)


# Plot 3: Directional change from first

mw_to_long_change_multi <- function(mw_list, model_name, params) {
  do.call(rbind, lapply(names(mw_list), function(pid) {
    df <- mw_list[[pid]]
    if (is.null(df)) return(NULL)

    keep <- intersect(params, names(df))
    if (length(keep) == 0) return(NULL)

    do.call(rbind, lapply(keep, function(p) {
      v <- df[[p]]
      idx0 <- which(is.finite(v))[1]  # first finite value
      if (length(idx0) == 0) return(NULL)
      v0 <- v[idx0]

      data.frame(
        Ppn = pid,
        model = model_name,
        start_trial = df$start_trial,
        parameter = p,
        change = v - v0,
        row.names = NULL
      )
    }))
  }))
}

df_long <- rbind(
  mw_to_long_change_multi(results$mw_linear, "Linear", params_to_check),
  mw_to_long_change_multi(results$mw_pt, "Prospect Theory", params_to_check)
)

# Summarize across participants per (model x start_trial x parameter)
df_sum <- aggregate(
  change ~ model + start_trial + parameter,
  data = df_long,
  FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
)

df_sum$mean <- df_sum$change[, "mean"]
df_sum$sd   <- df_sum$change[, "sd"]
df_sum$change <- NULL

df_sum$parameter <- factor(
  df_sum$parameter,
  levels = params_to_check,
  labels = c(expression(beta[PA]), expression(beta["NA"]),
             expression(gamma[PA]), expression(gamma["NA"]))
)

p3 <- ggplot(df_sum, aes(x = start_trial, y = mean, linetype = model)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.15) +
  geom_line(linewidth = 1) +
  facet_wrap(~ parameter, scales = "free_y", labeller = label_parsed) +
  labs(
    x = "Window start trial",
    y = "Change from first window",
    title = "Directional change in moving-window parameters across the task",
    subtitle = "Baseline-centered (0 = first window); line = mean; band = ±1 SD"
  ) +
  theme_word

print(p3)