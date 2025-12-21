# Objective function for VARMAX model estimation
objective_function_varmax <- function(parameters, data, restricted = FALSE) {
  
  # Cycle through participants 
  for (ppn in unique(data$Ppn)) {
    
    # Take the data by participant
    data_ppn <- subset(data, Ppn == ppn)
    
    # Load the data and split by participant
    data_ppn <- data[order(data$TrialNumber), ]
    y <- as.matrix(data_ppn[, c("PAscore", "NAscore")])  # d = 2
    x <- as.matrix(data_ppn[, c("Won", "TotalAmount")])  # m = 2
    Tn <- nrow(y)
    
    # Intercepts (2)
    alpha <- parameters[1:2]
    
    # Coefficient matrix (2x2 = 4 parameters)
    B_matrix <- matrix(parameters[3:6], nrow = 2, ncol = 2, byrow = TRUE)
    
    # Gamma matrix (2x2 = 4 parameters)
    Gamma_matrix <- matrix(parameters[7:10], nrow = 2, ncol = 2, byrow = TRUE)
    
    # Moving average matrix Phi (2x2 = 4 parameters)
    if (restricted) {
      # For restricted model: Phi = -Gamma
      Phi_matrix <- -Gamma_matrix
    } else {
      # For unrestricted model: estimate Phi separately
      Phi_matrix <- matrix(parameters[11:14], nrow = 2, ncol = 2, byrow = TRUE)
    }
    
    # For RESTRICTED model ( has fewer parameters):
    if (restricted) {
      cholesky_index <- 11  # Start at position 11
    } else {
      # For UNRESTRICTED model (has more parameters):
      cholesky_index <- 15  # Start at position 15
    }
    
    # Extract the 3 Cholesky parameters
    L <- matrix(0, 2, 2)
    L[lower.tri(L, diag = TRUE)] <- parameters[cholesky_index:(cholesky_index + 2)]
    
    # Build the covariance matrix
    Sigma <- L %*% t(L)
    
    # Stability conditions
    if (any(abs(eigen(Gamma_matrix)$values) >= 1)) {
      return(1e10)  # Return high penalty for unstable AR parameters
    }
    #Check for the MA parameters only if unrestricted
    if (!restricted && any(abs(eigen(Phi_matrix)$values) >= 1)) {
      return(1e10)  # Return high penalty for unstable MA parameters
    }
    
    # Delta 
    delta <- (diag(2) - Gamma_matrix) %*% alpha
    
    #Residuals and initial log-likelihood
    epsilon <- matrix(0, nrow = Tn, ncol = 2)
    nll <- 0
    
    # Start from t=2 since we need y_{t-1} and e_{t-1}
    for (t in 2:Tn) {
      # Predicted y at time t
      y_pred <- delta + B_matrix %*% x[t, ] + Gamma_matrix %*% y[t-1, ] + Phi_matrix %*% epsilon[t-1, ]
      
      # Residual at time t
      epsilon[t, ] <- y[t, ] - y_pred
      
      # Log-likelihood update
      nll <- nll - mvdnorm(epsilon[t, ], mean = c(0, 0), sigma = Sigma, log = TRUE)
    }
    
    total_nll <- total_nll + nll
  }
  return(total_nll)
}