library(DEoptim)

# Create data to estimate model to
data <- data.frame(x = rnorm(100))

# Objective function for the normal distribution
objective_function <- function(parameters) {
  # Get likelihoods for each datapoint
  L <- dnorm(
    data$x, 
    mean = parameters[1],
    sd = parameters[2]
  )
  
  # Transform to min-log-likelihood
  return(-sum(log(L)))
}

# Estimate normal model
results <- DEoptim(
  objective_function, 
  lower = c(-100, 10^(-5)),
  upper = c(100, 10^5),
  control = DEoptim.control(
    itermax = 100, 
    CR = 0.6,
    NP = 50,
    trace = FALSE
    # initialpop = 
  )
)

# Extract best parameters
results$optim$bestmem 

# Extract the MLL
results$optim$bestval
