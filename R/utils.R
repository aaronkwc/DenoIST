.denoist_cache_env <- new.env(parent = emptyenv())

denoist_get_offset_from_blocks <- function(idx, off_block_files, off_block_ranges) {
  block_idx <- which(idx >= off_block_ranges[, 1] & idx <= off_block_ranges[, 2])[1]
  if (length(block_idx) == 0 || is.na(block_idx)) {
    stop(sprintf("No offset block found for column index %d", idx))
  }
  blk <- readRDS(off_block_files[block_idx])
  local_idx <- idx - off_block_ranges[block_idx, 1] + 1L
  as.numeric(blk$data[, local_idx])
}


denoist_init_cluster <- function(requested_workers, prefer_fork = FALSE, verbose = FALSE) {
  if (requested_workers < 1L) {
    return(NULL)
  }

  types <- if (.Platform$OS.type == "windows") {
    c("PSOCK")
  } else if (prefer_fork) {
    c("FORK", "PSOCK")
  } else {
    c("PSOCK", "FORK")
  }

  for (tp in types) {
    workers <- as.integer(requested_workers)
    while (workers >= 1L) {
      clust <- tryCatch(
        parallel::makeCluster(workers, type = tp),
        error = function(e) {
          if (verbose) {
            message(sprintf("[debug] makeCluster(type=%s, workers=%d) failed: %s", tp, workers, e$message))
          }
          NULL
        }
      )

      if (!is.null(clust)) {
        return(list(cluster = clust, workers = workers, type = tp))
      }

      if (workers == 1L) {
        break
      }
      workers <- max(1L, workers %/% 2L)
    }
  }

  NULL
}


denoist_get_worker_rss <- function(clust) {
  rss_list <- parallel::clusterCall(clust, function() {
    pid <- Sys.getpid()
    rss_kb <- NA_real_
    status_path <- file.path("/proc", as.character(pid), "status")
    if (file.exists(status_path)) {
      status_lines <- readLines(status_path, warn = FALSE)
      rss_line <- status_lines[grepl("^VmRSS:", status_lines)]
      if (length(rss_line) > 0) {
        rss_kb <- suppressWarnings(as.numeric(sub(".*:\\s*([0-9]+)\\s*kB", "\\1", rss_line[1])))
      }
    }
    list(pid = pid, rss_kb = rss_kb)
  })

  pids <- vapply(rss_list, function(x) x$pid, numeric(1))
  rss_kb <- vapply(rss_list, function(x) x$rss_kb, numeric(1))
  list(pids = pids, rss_kb = rss_kb)
}


denoist_compact_param_result <- function(res, keep_posterior) {
  out <- list(
    lambda1 = res$lambda1,
    lambda2 = res$lambda2,
    pi = res$pi,
    log_lik = res$log_lik
  )
  if (keep_posterior) {
    out$posterior <- res$posterior
  }
  if (!is.null(res$error)) {
    out$error <- res$error
  }
  out
}


denoist_estimate_lambdas_moments <- function(x, s, pi) {
  keep <- s > 0
  if (!any(keep)) {
    return(c(lambda1 = NA_real_, lambda2 = NA_real_))
  }

  r <- x[keep] / s[keep]
  w <- s[keep]
  w_sum <- sum(w)
  if (!is.finite(w_sum) || w_sum <= 0) {
    return(c(lambda1 = NA_real_, lambda2 = NA_real_))
  }

  mean_r <- sum(w * r) / w_sum
  var_r <- sum(w * (r - mean_r)^2) / w_sum
  if (!is.finite(var_r) || var_r < 0) {
    var_r <- 0
  }

  if (!is.finite(pi) || pi <= 0 || pi >= 1) {
    return(c(lambda1 = mean_r, lambda2 = mean_r))
  }

  delta <- sqrt(var_r / (pi * (1 - pi)))
  lambda1 <- mean_r + delta
  lambda2 <- mean_r - delta
  if (!is.finite(lambda2) || lambda2 < 0) {
    lambda2 <- 0
  }

  c(lambda1 = lambda1, lambda2 = lambda2)
}


denoist_fixed_pi_posterior <- function(x, s, pi, lambda1, lambda2) {
  tau1 <- pi * dpois(x, s * lambda1)
  tau2 <- (1 - pi) * dpois(x, s * lambda2)
  gamma <- tau1 / (tau1 + tau2)
  gamma[!is.finite(gamma)] <- 0.5
  gamma
}
