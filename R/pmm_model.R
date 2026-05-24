#' @title Poisson mixture model solver
#' @param x A vector of counts.
#' @param s A vector of offsets.
#' @param max_iter Maximum number of iterations for the EM algorithm.
#' @param tol Tolerance for convergence.
#' @param n_inits Number of initial values for the mixing proportion inference. If input is a vector, directly use the vector has init values.
#' @param init_scheme Initialization scheme for the EM algorithm. Options are 'random' for random initialization and 'kmeanspp' for k-means++ initialization.
#' @param seed Random seed used for k-means++ initialization.
#' @param posterior_cutoff Cutoff for posterior probabilities to assign memberships.
#' @param return_posterior Logical, if TRUE (default), return full posterior probabilities. Set to FALSE to reduce memory usage.
#' @param use_cpp Logical, if TRUE, use the Rcpp EM loop for faster in-place computations.
#' @param membership_output Format of membership return value. Options are 'full' (default) to return a full 0/1 vector, and 'sparse' to return only zero-membership indices.
#' @param pad_length Logical, if TRUE, pad outputs to the original input length. Currently kept for backward compatibility.
#' @param verbose Logical, if TRUE, print progress messages.
#' @return{
#' A list containing the following elements:
#' \item{memberships}{A vector indicating the membership of each observation (0 or 1), returned when `membership_output = "full"`.}
#' \item{membership_zero_idx}{Integer indices where membership is 0, returned when `membership_output = "sparse"`.}
#' \item{n}{Original input length, returned when `membership_output = "sparse"`.}
#' \item{posterior}{A vector of posterior probabilities for each observation, returned when `return_posterior = TRUE`.}
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
#' @export
solve_poisson_mixture <- function(x, s,
                                  max_iter = 5000,
                                  tol = 1e-6,
                                  n_inits = 10,
                                  init_scheme = 'random',
                                  seed = 42, # only used for kmeans++ init
                                  posterior_cutoff = 0.6,
                                  return_posterior = TRUE,
                                  use_cpp = TRUE,
                                  membership_output = c("full", "sparse"),
                                  pad_length = FALSE,
                                  verbose = FALSE) {

  membership_output <- match.arg(membership_output)

  n <- length(x)

  format_membership <- function(full_memberships) {
    full_memberships <- as.integer(full_memberships)
    if (membership_output == "sparse") {
      return(list(membership_zero_idx = which(full_memberships == 0L),
                  n = length(full_memberships)))
    }
    list(memberships = full_memberships)
  }

  build_output <- function(full_memberships,
                           lambda1,
                           lambda2,
                           pi,
                           log_lik,
                           full_posterior = NULL) {
    out <- format_membership(full_memberships)
    if (return_posterior) {
      out$posterior <- full_posterior
    }
    out$lambda1 <- lambda1
    out$lambda2 <- lambda2
    out$pi <- pi
    out$log_lik <- log_lik
    out
  }

  # Degenerate edge case: empty input vector.
  if(n == 0L){
    return(build_output(full_memberships = integer(0),
                        full_posterior = numeric(0),
                        lambda1 = NA_real_,
                        lambda2 = NA_real_,
                        pi = NA_real_,
                        log_lik = NA_real_))
  }

  if(length(n_inits) == 1) {
    # if n_inits is a single number greater than or equal to 1, generate that many random initial values for pi
    if(n_inits >= 1){pi_inits <- runif(n_inits, min = 0, max = 0.5)}else if(n_inits < 1){
      # directly use the 1 init
      pi_inits <- n_inits
    }
  }else{
    pi_inits <- n_inits
  }

  # check that pi_inits are all >0 and <=1
  if(any(pi_inits <= 0) || any(pi_inits > 1)) {
    stop("All initial values for pi must be in the range (0, 1]")
  }

  # Store indices of non-zero s
  non_zero_indices <- which(s > 0)

  # Degenerate edge case: no valid offsets remain after filtering.
  if(length(non_zero_indices) == 0L){
    return(build_output(full_memberships = rep.int(1L, n),
                        full_posterior = rep.int(1, n),
                        lambda1 = NA_real_,
                        lambda2 = NA_real_,
                        pi = NA_real_,
                        log_lik = NA_real_))
  }

  # Filter to entries with non-zero offsets for model fitting.
  x_filtered <- x[non_zero_indices]
  s_filtered <- s[non_zero_indices]

  # Degenerate edge case: all observed counts are zero.
  if(all(x_filtered == 0)){
    return(build_output(full_memberships = rep.int(1L, n),
                        full_posterior = rep.int(1, n),
                        lambda1 = 0,
                        lambda2 = 0,
                        pi = 1,
                        log_lik = 0))
  }

  best_result <- NULL
  best_log_lik <- -Inf

  if(init_scheme == 'kmeanspp'){
    if(verbose){
      message("Using k-means++ for initialization...")
    }
    init_res <- kmeanspp_poisson_init(x_filtered, s_filtered, K = 2, nstart = 1, seed = seed)
    pi_inits <- min(init_res$pi)
    lambda1_init <- max(init_res$lambda)
    lambda2_init <- min(init_res$lambda)
  }else if(init_scheme == 'random'){
    if(verbose){
      message("Using random initialization...")
    }
    lambda1_init <- NA
    lambda2_init <- NA
  }else{
    stop("Invalid init_scheme. Choose 'random' or 'kmeanspp'.")
  }

  if (!use_cpp) {
    tau1 <- numeric(length(x_filtered))
    tau2 <- numeric(length(x_filtered))
    gamma <- numeric(length(x_filtered))
  }

  for (pi_init in pi_inits) {
    # Initialize parameters
    if(is.na(lambda1_init)){lambda1 <- mean(x_filtered) / mean(s_filtered)}else{lambda1 <- lambda1_init}
    if(is.na(lambda2_init)){lambda2 <- mean(x_filtered) / (2 * mean(s_filtered))}else{lambda2 <- lambda2_init}

    pi <- pi_init

    if (verbose) {
      message("Initial parameters for pi =", pi_init, ":\n")
      message("lambda1:", lambda1, "lambda2:", lambda2, "pi:", pi, "\n")
    }

    if (use_cpp) {
      if (verbose) {
        message("Running C++ EM loop...")
      }
      cpp_res <- poisson_mixture_em_cpp(x_filtered, s_filtered, lambda1, lambda2, pi, max_iter, tol)
      lambda1 <- cpp_res$lambda1
      lambda2 <- cpp_res$lambda2
      pi <- cpp_res$pi
      log_lik <- cpp_res$log_lik
      gamma <- cpp_res$gamma
    } else {
      log_likelihood <- function(x, s, lambda1, lambda2, pi) {
        sum(log(pi * dpois(x, s * lambda1) + (1 - pi) * dpois(x, s * lambda2)))
      }

      log_lik <- log_likelihood(x_filtered, s_filtered, lambda1, lambda2, pi)

      for (iter in seq_len(max_iter)) {
        # E-step: calculate responsibilities
        tau1[] <- pi * dpois(x_filtered, s_filtered * lambda1)
        tau2[] <- (1 - pi) * dpois(x_filtered, s_filtered * lambda2)
        gamma[] <- tau1 / (tau1 + tau2)
        # Degenerate edge case: both components underflow to zero probability.
        gamma[!is.finite(gamma)] <- 0.5

        # M-step: update parameters
        gamma_s <- sum(gamma * s_filtered)
        #one_minus_gamma <- 1 - gamma
        one_minus_gamma_s <- sum(1-gamma * s_filtered)

        # Degenerate edge case: M-step denominator collapses to zero.
        if(gamma_s <= 0 || one_minus_gamma_s <= 0){
          break
        }

        lambda1 <- sum(gamma * x_filtered) / gamma_s
        lambda2 <- sum((1-gamma) * x_filtered) / one_minus_gamma_s
        pi <- mean(gamma)

        if (verbose) {
          message("Iteration", iter, "parameters:\n")
          message("lambda1:", lambda1, "lambda2:", lambda2, "pi:", pi, "\n")
        }

        # Check for convergence
        new_log_lik <- log_likelihood(x_filtered, s_filtered, lambda1, lambda2, pi)
        if (!is.finite(new_log_lik) || abs(new_log_lik - log_lik) < tol) {
          if (verbose) {
            message("Converged after", iter, "iterations\n")
            message("Final log-likelihood:", log_lik, "\n")
          }
          break
        }

        log_lik <- new_log_lik
      }
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
                gamma = rep(1, length(x_filtered)))
      }
    }
  }

  # Assign memberships
  memberships <- ifelse(best_result$gamma >= posterior_cutoff, 1, 0)

  # TODO: if memberships are all 0, set to 1
  if (all(memberships == 0)) {
    memberships <- rep(1, length(memberships))
  }

  # Pad the results to match the original input length
  full_memberships <- rep.int(1L, n)
  full_memberships[non_zero_indices] <- memberships

  full_posterior <- NULL
  if(return_posterior){
    full_posterior <- rep.int(1, n)
    full_posterior[non_zero_indices] <- best_result$gamma
  }

  return(build_output(full_memberships = full_memberships,
                      full_posterior = full_posterior,
                      lambda1 = best_result$lambda1,
                      lambda2 = best_result$lambda2,
                      pi = best_result$pi,
                      log_lik = best_result$log_lik))
}


kmeanspp_poisson_init <- function(x, offset, K = 2, nstart = 1, seed = NULL) {
  if(!is.null(seed)) set.seed(seed)

  # Step 0: compute rate-normalized counts
  r <- x / offset

  # Step 1: choose the first center randomly
  centers <- numeric(K)
  centers[1] <- sample(r, 1)

  # Step 2: choose remaining centers using probability proportional to squared distance
  for(k in 2:K) {
    # squared distance to nearest existing center
    dist2 <- sapply(r, function(v) min((v - centers[1:(k-1)])^2))
    prob <- dist2 / sum(dist2)
    centers[k] <- sample(r, 1, prob = prob)
  }

  # Step 3: assign each point to nearest center
  assign <- sapply(r, function(v) which.min((v - centers)^2))

  # Step 4: compute initial lambdas and pis
  lambda_init <- numeric(K)
  pi_init <- numeric(K)

  for(k in 1:K) {
    if(sum(assign == k) == 0) {
      # fallback if no points assigned
      lambda_init[k] <- centers[k]
      pi_init[k] <- 1/K
    } else {
      lambda_init[k] <- sum(x[assign == k]) / sum(offset[assign == k])
      pi_init[k] <- sum(assign == k) / length(x)
    }
  }

  # Step 5: order lambdas
  ord <- order(lambda_init)
  lambda_init <- lambda_init[ord]
  pi_init <- pi_init[ord]

  list(pi = pi_init, lambda = lambda_init)

  }
