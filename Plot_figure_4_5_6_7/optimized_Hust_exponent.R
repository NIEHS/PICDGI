
    ################# Optimize Hurst exponent H ---------------
    ################# 

    library(Matrix)
    library(ggplot2)
    
    # Function to compute covariance matrix for given H
    compute_covariance <- function(H, N, nu.square, coef) {
      par = c(H, nu.square, coef)
      cor <- sapply(1:N, function(j) 0.5 * (((abs(j+1))^(2.*par[1])) - 
                                              2*((abs(j))^(2.*par[1])) + 
                                              ((abs(j-1))^(2.*par[1]))))
      cor.matrix <- toeplitz(cor)
      return(cor.matrix)
    }
    
    # Parameters
    N <- 100
    nu.square <- 1.78
    coef <- 1.15
    
    # Generate covariance matrices for different H values
    H_values <- seq(0.1, 0.9, by=0.2)
    cov_matrices <- lapply(H_values, function(H) compute_covariance(H, N, nu.square, coef))
    
    # Plot Covariance Matrices for Different H
    par(mfrow=c(2,3))
    for (i in 1:length(H_values)) {
      image(cov_matrices[[i]], main = paste("H =", H_values[i]), col=heat.colors(10))
    }
    
    # Example observed covariance matrix (simulated real data)
    observed_cov <- compute_covariance(0.6, N, nu.square, coef)
    
    # Function to compute error between estimated and observed covariance matrices
    compute_error <- function(H, N, nu.square, coef, observed_cov) {
      estimated_cov <- compute_covariance(H, N, nu.square, coef)
      sum((estimated_cov - observed_cov)^2)  # Sum of squared differences
    }
    
    # Optimize H by minimizing the error
    H_range <- seq(0.1, 0.9, by=0.01)
    errors <- sapply(H_range, function(H) compute_error(H, N, nu.square, coef, observed_cov))
    optimal_H <- H_range[which.min(errors)]
    
    # Plot the optimization process
    error_data <- data.frame(H = H_range, Error = errors)
    
    ggplot(error_data, aes(x = H, y = Error)) +
      geom_line(color = "blue", size = 1) +
      geom_point(aes(x = optimal_H, y = min(errors)), color = "red", size = 10) +
      ggtitle("Optimization of H") +
      xlab("H Value") +
      ylab("Error") +
      theme_minimal()
    
    # Print the optimal H
    print(paste("Optimal H:", optimal_H))

    
    
    
    
        
    
    
    
    
  
  
  
  
  
  