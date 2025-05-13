#' @title Poisson mixture model solver
#' @param x A vector of counts.
#' @param s A vector of offsets.
#' @param max_iter Maximum number of iterations for the EM algorithm.
#' @param tol Tolerance for convergence.
#' @param pi_inits Initial values for the mixing proportion.
#' @param verbose Logical, if TRUE, print progress messages.
#' @return{
#' A list containing the following elements:
#' \item{memberships}{A vector indicating the membership of each observation (0 or 1).}
#' \item{posterior}{A vector of posterior probabilities for each observation.}
#' \item{lambda1}{Estimated parameter for the first component.}
#' \item{lambda2}{Estimated parameter for the second component.}
#' \item{pi}{Estimated mixing proportion.}
#' \item{log_lik}{Log-likelihood of the fitted model.}
#' }
#' @description This function implements a 2-component Poisson mixture model
#' using the EM algorithm. It estimates the parameters of the model and assigns
#' memberships to each observation based on the posterior probabilities.
#' @details The function takes a vector of counts and a vector of offsets as input.
#' It uses the EM algorithm to iteratively update the parameters of the model
#' until convergence is reached or the maximum number of iterations is exceeded.
#' The function also allows for multiple initialisations of the mixing proportion
#' to find the best solution.
#' @examples
#' x <- rpois(100, lambda = 5)
#' s <- runif(100, min = 0, max = 1)
#' result <- solve_poisson_mixture(x, s)
#' print(result)
#' @import stats
#' @export
solve_poisson_mixture <- function(x, s,
                                  max_iter = 5000,
                                  tol = 1e-6,
                                  pi_inits = runif(10, min = 0, max = 0.5),
                                  verbose = FALSE) {

  n <- length(x)

  # Store indices of non-zero s
  non_zero_indices <- which(s > 0)

  # Remove entries with s = 0
  x <- x[non_zero_indices]
  s <- s[non_zero_indices]

  best_result <- NULL
  best_log_lik <- -Inf

  for (pi_init in pi_inits) {
    # Initialize parameters
    lambda1 <- mean(x) / mean(s)
    lambda2 <- mean(x) / (2 * mean(s))
    pi <- pi_init

    if (verbose) {
      cat("Initial parameters for pi =", pi_init, ":\n")
      cat("lambda1:", lambda1, "lambda2:", lambda2, "pi:", pi, "\n")
    }

    log_likelihood <- function(x, s, lambda1, lambda2, pi) {
      sum(log(pi * dpois(x, s * lambda1) + (1 - pi) * dpois(x, s * lambda2)))
    }

    log_lik <- log_likelihood(x, s, lambda1, lambda2, pi)

    for (iter in 1:max_iter) {
      # E-step: calculate responsibilities
      tau1 <- pi * dpois(x, s * lambda1)
      tau2 <- (1 - pi) * dpois(x, s * lambda2)
      gamma <- tau1 / (tau1 + tau2)

      # M-step: update parameters
      lambda1 <- sum(gamma * x) / sum(gamma * s)
      lambda2 <- sum((1 - gamma) * x) / sum((1 - gamma) * s)
      pi <- mean(gamma)

      if (verbose) {
        cat("Iteration", iter, "parameters:\n")
        cat("lambda1:", lambda1, "lambda2:", lambda2, "pi:", pi, "\n")
      }

      # Check for convergence
      new_log_lik <- log_likelihood(x, s, lambda1, lambda2, pi)
      if (!is.finite(new_log_lik) || abs(new_log_lik - log_lik) < tol) {
        if (verbose) {
          cat("Converged after", iter, "iterations\n")
          cat("Final log-likelihood:", log_lik, "\n")
        }
        break
      }

      log_lik <- new_log_lik
    }

    if (log_lik > best_log_lik) {
      best_log_lik <- log_lik
      if(abs(lambda1 - lambda2) > 1e-2) {
        # Store the best parameters
        best_result <- list(lambda1 = lambda1,
                            lambda2 = lambda2,
                            pi = pi,
                            log_lik = log_lik,
                            gamma = gamma)
      }else{
        # If model collapse occurs, keep everything
        best_result <- list(lambda1 = lambda1,
                            lambda2 = lambda2,
                            pi = pi,
                            log_lik = log_lik,
                            gamma = rep(1, length(x)))
      }
    }
  }

  # Assign memberships
  memberships <- ifelse(best_result$gamma > 0.6, 1, 0)

  # TODO: if memberships are all 0, set to 1
  if (all(memberships == 0)) {
    memberships <- rep(1, length(memberships))
  }

  # Pad the results to match the original input length
  full_memberships <- rep(1, n)
  full_memberships[non_zero_indices] <- memberships

  full_posterior <- rep(1, n)
  full_posterior[non_zero_indices] <- best_result$gamma

  return(list(memberships = full_memberships,
              posterior = full_posterior,
              lambda1 = best_result$lambda1,
              lambda2 = best_result$lambda2,
              pi = best_result$pi,
              log_lik = best_result$log_lik))
}

# Function to apply solve_poisson_mixture to a single column
apply_poisson_mixture_single <- function(index, c_matrix, s_matrix) {
  test_c <- c_matrix[, index]
  test_s <- s_matrix[, index]
  result <- tryCatch({
    solve_poisson_mixture(test_c, test_s)
  }, error = function(e) {
    list(memberships = rep(1, length(test_c)),
         posterior = rep(1, length(test_c)),
         lambda1 = NA, lambda2 = NA, pi = NA)
  })
  return(result)
}
