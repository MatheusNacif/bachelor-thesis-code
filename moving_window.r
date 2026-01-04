
# Moving-window test for T = 152 trials (ppn 1 only)


library(DEoptim)
library(ggplot2)


data <- read.csv("Data_CIAC.csv")

# Keep only first participant and first 152 trials
ppn_id <- data$Ppn[1]
data_ppn <- data %>%
  filter(Ppn == ppn_id) %>%
  arrange(TrialNumber)

if (nrow(data_ppn) < 152) stop("Participant has fewer than 152 trials.")
data_ppn <- data_ppn[1:152, ]

Tn <- nrow(data_ppn)
cat("Using Ppn =", ppn_id, "with", Tn, "trials\n")

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

# Bounds
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

p_vec = length(lower)

# Window sizes to test + sliding step
window_sizes <- c(15, 20, 25, 30, 35, 40, 50, 60, 70, 75, 80, 90, 100)
step_size <- 1

# DEoptim control
de_ctrl <- DEoptim.control(
  NP = p_vec * 10,          
  itermax = 80,     
  trace = FALSE
)

# Fit model for one window
fit_one_window <- function(data_window) {
  fit <- DEoptim(
    fn = objective_exp_discount_nll,
    lower = lower,
    upper = upper,
    data_ppn = data_window,
    control = de_ctrl
  )
  list(
    par = fit$optim$bestmem,
    obj = fit$optim$bestval
  )
}

# Run moving-window fits for a given window size 
run_for_window_size <- function(w) {
  starts <- seq(1, Tn - w + 1, by = step_size)

  pars_mat <- matrix(NA_real_, nrow = length(starts), ncol = length(lower))
  obj_vec  <- rep(NA_real_, length(starts))

  for (i in seq_along(starts)) {
    s <- starts[i]
    e <- s + w - 1
    data_window <- data_ppn[s:e, ]

    out <- fit_one_window(data_window)
    pars_mat[i, ] <- out$par
    obj_vec[i] <- out$obj

    cat(sprintf("w=%3d  window %3d/%3d  trials [%3d,%3d]  obj=%.3f\n",
                w, i, length(starts), s, e, out$obj))
  }

  colnames(pars_mat) <- paste0("p", seq_len(ncol(pars_mat)))

  list(
    w = w,
    starts = starts,
    pars = pars_mat,
    obj = obj_vec
  )
}

# Instability metric across windows
# Here: coefficient of variation-like measure per parameter:
#   instability_k = sd(par_k) / (|mean(par_k)| + eps)
# Then average across parameters.
instability_metric <- function(pars_mat, eps = 1e-6) {
  m <- colMeans(pars_mat, na.rm = TRUE)
  s <- apply(pars_mat, 2, sd, na.rm = TRUE)
  mean(s / (abs(m) + eps), na.rm = TRUE)
}

# Run all window sizes
all_results <- lapply(window_sizes, run_for_window_size)

# Summarize per window size
summary_df <- bind_rows(lapply(all_results, function(res) {
  data.frame(
    window_size = res$w,
    n_windows = length(res$obj),
    mean_obj = mean(res$obj, na.rm = TRUE),
    # normalize by number of trials in window to compare across w
    mean_obj_per_trial = mean(res$obj / res$w, na.rm = TRUE),
    instability = instability_metric(res$pars)
  )
}))

print(summary_df)

# Pick a "sweet spot" with a combined score 
# Score = mean_obj_per_trial + lambda * instability
# Increase lambda if you care more about stability than fit.
lambda <- 0.25
summary_df <- summary_df %>%
  mutate(score = mean_obj_per_trial + lambda * instability) %>%
  arrange(score)

cat("\nTop 5 window sizes by combined score (lower is better):\n")
print(head(summary_df, 5))

sweetspot <- summary_df$window_size[1]
cat("\n Sweet spot window size:", sweetspot, "\n")

# Plots
p1 <- ggplot(summary_df, aes(x = window_size, y = mean_obj_per_trial)) +
  geom_line() + geom_point() +
  labs(title = "Average objective per trial vs window size",
       x = "Window size (trials)", y = "Mean objective / trial")

p2 <- ggplot(summary_df, aes(x = window_size, y = instability)) +
  geom_line() + geom_point() +
  labs(title = "Parameter instability vs window size",
       x = "Window size (trials)", y = "Instability (avg sd/|mean|)")

p3 <- ggplot(summary_df, aes(x = window_size, y = score)) +
  geom_line() + geom_point() +
  labs(title = paste0("Combined score vs window size (lambda=", lambda, ")"),
       x = "Window size (trials)", y = "Score (lower is better)")

print(p1); print(p2); print(p3)
