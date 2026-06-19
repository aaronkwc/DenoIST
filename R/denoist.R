utils::globalVariables(c("feature_name", "hexbin_id", "count"))

#' @importFrom pbapply pblapply
#' @importFrom stats dpois runif
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment assay
#' @importFrom parallel makeCluster clusterExport parLapply stopCluster clusterEvalQ clusterCall
#' @importFrom HDF5Array writeHDF5Array HDF5Array
#' @importFrom DelayedArray path
#' @useDynLib DenoIST, .registration = TRUE
#' @export
#' @title DenoIST
#' @description
#' DenoIST (Denoising Image-based Spatial Transcriptomics) is a method for
#' identifying and removing contamination artefacts in
#' image-based single-cell transcriptomics (IST) data.
#' It uses a transposed Poisson mixture model to identify contamination.
#' @param mat A matrix of counts (genes x cells), or a SpatialExperiment object.
#' @param coords A data frame of coordinates (n_cells x 2). Only used if not using a SpatialExperiment object as input.
#' @param tx A data frame of transcript with x, y and qv columns.
#' @param tx_x Column name for the x coordinates in the transcripts dataframe. Default is 'x'.
#' @param tx_y Column name for the y coordinates in the transcripts dataframe. Default is 'y'.
#' @param feature_label Column name for the gene of each transcript in the transcripts dataframe. Default is 'gene'.
#' @param distance The maximum distance to consider for local background estimation.
#' @param nbins The number of bins to use for hexagonal binning.
#' @param posterior_cutoff The cutoff for posterior probability to determine contamination.
#' @param n_inits The number of initialisations for the mixing proportion in the Poisson mixture model. If input is a vector, directly use the vector of values as init values.
#' @param cl The number of cores to use for parallel processing.
#' @param out_dir The output directory to save the results.
#' @param neighbour_mode The method to use for finding neighbors. Options are 'fast' (default) and 'original'. 'fast' uses a faster implementation with dbscan, while 'original' uses the original implementation.
#' @param init_scheme The method to use for initializing the mixing proportion in the Poisson mixture model. Options are 'random' (default) and 'kmeans++'. 'random' generates random initial values, while 'kmeans++' uses a k-means++ inspired approach to generate initial values.
#' @param on_disk Logical, if TRUE, use on-disk processing with HDF5Array to handle large datasets that may not fit in memory. If FALSE (default), process in-memory.
#' @param offset_on_disk Logical, if TRUE, compute neighbour offsets in on-disk block mode for the fast neighbour function. If FALSE (default), keep the original in-memory off_mat behavior.
#' @param block_size The number of cells to process in each block when using parallel processing. Adjust this based on available memory and number of cores.
#' @param return_params Logical, if TRUE (default), return per-cell model parameters in the `params` output slot. Set to FALSE to reduce memory usage.
#' @param keep_posterior Logical, if TRUE (default), include posterior vectors in `params`. Set to FALSE to reduce memory usage when processing large datasets.
#' @param sparse_membership Logical, if TRUE, request sparse per-cell membership output from the PMM solver during computation and return memberships in an index-based sparse representation (zero indices per column), avoiding dense memberships matrix allocation.
#' @param verbose Logical, if TRUE, print progress messages.
#' @return A list containing the following elements:
#' \item{memberships}{A matrix of memberships for each gene in each cell.}
#' \item{adjusted_counts}{A matrix of adjusted counts for each gene in each cell.}
#' \item{params}{A list of parameters for each gene.}
#' @details
#' The function calculates local background using hexagonal binning and applies
#' a Poisson mixture model to identify contamination.
#' It returns a matrix of memberships and adjusted counts for each gene in each cell.
#' @examples
#' # Load example data
#' set.seed(42)
#' mat <- matrix(rpois(1000, lambda = 10), nrow = 10, ncol = 100)
#' rownames(mat) <- paste0("gene", 1:10)
#' coords <- data.frame(x = rnorm(100), y = rnorm(100))
#' tx <- data.frame(x = c(rnorm(500), rnorm(500, 3)),
#'                  y = c(rnorm(500), rnorm(500, 3)),
#'                  qv = rep(30, 1000), gene = paste0('gene', 1:10))
#' # Run DenoIST
#' result <- denoist(mat, tx, coords, distance = 1, nbins = 50, cl = 1,
#'                   out_dir = NULL, verbose = TRUE)
#' # Check results
#' print(result$memberships[1:5, 1:5])
#' print(result$adjusted_counts[1:5, 1:5])
#' print(result$params[[1]])
denoist <- function(mat, tx, coords = NULL,
                    tx_x = "x",
                    tx_y = "y",
                    feature_label = "gene",
                    distance = 50,
                    nbins = 200,
                    posterior_cutoff = 0.6,
                    n_inits = 10,
                    cl = 1,
                    out_dir = NULL,
                    neighbour_mode = 'fast',
                    init_scheme = 'random',
                    on_disk = FALSE,
                    offset_on_disk = FALSE,
                    block_size = 5000,
                    return_params = TRUE,
                    keep_posterior = TRUE,
                    sparse_membership = FALSE,
                    verbose = FALSE){

  # -------------------------------
  # Input handling
  # -------------------------------
  if(inherits(mat, "SpatialExperiment")){
    coords <- spatialCoords(mat)
    mat <- assay(mat)
    mat <- mat[!grepl("NegControl|BLANK|Unassigned", rownames(mat)),]
  } else if(is.null(coords)){
    stop("coords must be provided")
  }

  # -------------------------------
  # Compute neighbour offset
  # -------------------------------
  if(verbose) message("Calculating neighbour offset...")

  if(neighbour_mode == 'fast'){
    off_mat <- local_offset_distance_with_background_fast(
      mat = mat, coords = coords, tx = tx,
      tx_x = tx_x, tx_y = tx_y,
      feature_label = feature_label,
      distance = distance, nbins = nbins,
      cl = cl, verbose = verbose,
      on_disk = offset_on_disk,
      block_size = block_size
    )
  } else if(neighbour_mode == 'original'){
    off_mat <- local_offset_distance_with_background(
      mat = mat, coords = coords, tx = tx,
      tx_x = tx_x, tx_y = tx_y,
      feature_label = feature_label,
      distance = distance, nbins = nbins,
      cl = cl, verbose = verbose
    )
  } else {
    stop("Invalid neighbour_mode. Choose 'fast' or 'original'.")
  }


  # free memory
  rm(tx)
  gc()

  is_off_blocks <- is.list(off_mat) && !is.null(off_mat$mode) && identical(off_mat$mode, "on_disk_blocks")

  get_offset_col <- if (is_off_blocks) {
    off_block_files <- off_mat$files
    off_block_ranges <- off_mat$block_col_ranges
    function(idx) denoist_get_offset_from_blocks(idx, off_block_files, off_block_ranges)
  } else {
    function(idx) off_mat[, idx]
  }

  pmm_return_posterior <- isTRUE(return_params) && isTRUE(keep_posterior)
  pmm_membership_output <- if (isTRUE(sparse_membership)) "sparse" else "full"
  use_sparse_memberships_out <- isTRUE(sparse_membership)

  get_count_col <- function(idx) mat[, idx]

  # Function to parallelize over
  apply_poisson_mixture_single <- function(get_count_col,
                                           get_offset_col,
                                           posterior_cutoff,
                                           n_inits,
                                           init_scheme,
                                           pmm_return_posterior,
                                           pmm_membership_output) {
    function(idx) { # friendship ended with OOP, FP is now my best friend
      c_vec <- get_count_col(idx)
      s_vec <- get_offset_col(idx)

      result <- tryCatch({
        out <- solve_poisson_mixture(c_vec, s_vec,
                                     posterior_cutoff = posterior_cutoff,
                                     n_inits = n_inits,
                                     init_scheme = init_scheme,
                                     return_posterior = pmm_return_posterior,
                                     membership_output = pmm_membership_output,
                                     use_cpp = TRUE,
                                     max_iter = 500)
        return(out)
      }, error = function(e) {
        out <- list(lambda1 = NA, lambda2 = NA, pi = NA, log_lik = NA, error = e$message)
        if(pmm_membership_output == "sparse"){
          out$membership_zero_idx <- integer(0)
          out$n <- length(c_vec)
        } else {
          out$memberships <- rep.int(1L, length(c_vec))
        }
        if(pmm_return_posterior){
          out$posterior <- rep.int(1, length(c_vec))
        }
        return(out)
      })
      return(result)
    }
  }

  memberships_matrix <- if(!use_sparse_memberships_out){
    matrix(1,
           nrow = nrow(mat),
           ncol = ncol(mat),
           dimnames = list(rownames(mat), colnames(mat)))
  } else {
    NULL
  }
  memberships_zero_idx_by_col <- if(use_sparse_memberships_out){
    vector("list", ncol(mat))
  } else {
    NULL
  }
  params_out <- if(return_params) vector("list", ncol(mat)) else NULL

  store_results <- function(indices, res_list) {
    for(k in seq_along(indices)){
      idx <- indices[k]
      res <- res_list[[k]]
      if(use_sparse_memberships_out){
        zero_idx <- if(!is.null(res$membership_zero_idx)){
          as.integer(res$membership_zero_idx)
        } else if(!is.null(res$memberships)){
          which(res$memberships == 0L)
        } else {
          stop("Invalid PMM result: expected memberships or membership_zero_idx")
        }
        memberships_zero_idx_by_col[[idx]] <<- zero_idx
      } else {
        if(!is.null(res$memberships)){
          memberships_matrix[, idx] <<- res$memberships
        } else if(!is.null(res$membership_zero_idx)){
          if(length(res$membership_zero_idx) > 0L){
            memberships_matrix[res$membership_zero_idx, idx] <<- 0L
          }
        } else {
          stop("Invalid PMM result: expected memberships or membership_zero_idx")
        }
      }
      if(return_params){
        params_out[[idx]] <<- denoist_compact_param_result(res, keep_posterior)
      }
    }
  }

  # -------------------------------
  # Parallel / sequential execution
  # -------------------------------
  if(cl > 1){
    if(on_disk){
      if(verbose) message("Using on-disk processing with HDF5Array...")
      # -------------------------------
      # Write matrices to HDF5
      # -------------------------------
      if(verbose) message("Writing matrices to HDF5...")

      mat_tmp_path <- tempfile("mat_", fileext = ".h5")

      mat_h5 <- writeHDF5Array(mat, filepath = mat_tmp_path,
                               name = "counts",
                               chunkdim = c(nrow(mat), 1))

      use_offset_blocks <- is_off_blocks
      if(!use_offset_blocks){
        off_tmp_path <- tempfile("off_", fileext = ".h5")
        off_h5 <- writeHDF5Array(off_mat, filepath = off_tmp_path,
                                 name = "offsets",
                                 chunkdim = c(nrow(off_mat), 1))
      } else {
        off_block_files <- off_mat$files
        off_block_ranges <- off_mat$block_col_ranges
      }

      # free RAM copies not needed beyond this point for offset
      rm(off_mat)
      gc()

      h5_mat_path <- path(mat_h5)
      if(!use_offset_blocks){
        h5_off_path <- path(off_h5)
      }

      n_cells <- ncol(mat_h5)

        # Register cleanup immediately
        on.exit({
          if (file.exists(h5_mat_path)) unlink(h5_mat_path, force = TRUE)
          if (!use_offset_blocks && file.exists(h5_off_path)) unlink(h5_off_path, force = TRUE)
        }, add = TRUE)

        cluster_info <- denoist_init_cluster(cl - 1L, prefer_fork = FALSE, verbose = verbose)
        clust <- NULL
        if(is.null(cluster_info)){
          warning("Cluster setup failed for all worker counts/types; falling back to sequential PMM execution.")
        } else {
          clust <- cluster_info$cluster
          if(verbose){
            message(sprintf("Using %d worker(s) with %s cluster.", cluster_info$workers, cluster_info$type))
          }
        }

        if(verbose) message("Exporting variables to workers...")

        # -------------------------------
        # Worker function
        # -------------------------------
        worker_fun <- function(idx) {
          c_vec <- as.numeric(mat_h5[, idx])
          if(use_offset_blocks){
            block_idx <- findInterval(idx, off_block_starts)
            if(block_idx < 1L || block_idx > length(off_block_files) || idx > off_block_ends[block_idx]){
              stop(sprintf("No offset block found for column index %d", idx))
            }
            cache <- get0(".denoist_offset_cache", envir = .denoist_cache_env, inherits = FALSE)
            if(is.null(cache) || cache$block_idx != block_idx){
              blk <- readRDS(off_block_files[block_idx])
              cache <- list(block_idx = block_idx, data = blk$data)
              assign(".denoist_offset_cache", cache, envir = .denoist_cache_env)
            }
            local_idx <- idx - off_block_starts[block_idx] + 1L
            s_vec <- as.numeric(cache$data[, local_idx])
          } else {
            s_vec <- as.numeric(off_h5[, idx])
          }

          tryCatch({
            solve_poisson_mixture(
              c_vec, s_vec,
              posterior_cutoff = posterior_cutoff,
              n_inits = n_inits,
              init_scheme = init_scheme,
              return_posterior = pmm_return_posterior,
              membership_output = pmm_membership_output
            )
          }, error = function(e) {
            out <- list(
              lambda1 = NA, lambda2 = NA, pi = NA,
              log_lik = NA,
              error = e$message
            )
            if(pmm_membership_output == "sparse"){
              out$membership_zero_idx <- integer(0)
              out$n <- length(c_vec)
            } else {
              out$memberships <- rep.int(1L, length(c_vec))
            }
            if(pmm_return_posterior){
              out$posterior <- rep.int(1, length(c_vec))
            }
            out
          })
        }


        if(use_offset_blocks){
          off_block_starts <- off_block_ranges[, 1]
          off_block_ends <- off_block_ranges[, 2]
        }

        if(!is.null(clust)){
          clusterExport(clust, c(
            "worker_fun",
            "solve_poisson_mixture",
            "posterior_cutoff",
            "n_inits",
            "init_scheme",
            "pmm_return_posterior",
            "pmm_membership_output",
            "h5_mat_path",
            "use_offset_blocks"
          ), envir = environment())

          if(use_offset_blocks){
            clusterExport(clust, c("off_block_files", "off_block_starts", "off_block_ends"), envir = environment())
          } else {
            clusterExport(clust, c("h5_off_path"), envir = environment())
          }

          # Initialize workers ONCE
          clusterEvalQ(clust, {
            mat_h5 <- HDF5Array::HDF5Array(h5_mat_path, "counts")
            if(!use_offset_blocks){
              off_h5 <- HDF5Array::HDF5Array(h5_off_path, "offsets")
            } else {
              assign(".denoist_offset_cache", NULL, envir = .denoist_cache_env)
            }
            NULL
          })
        } else if(use_offset_blocks){
          assign(".denoist_offset_cache", NULL, envir = .denoist_cache_env)
        }

        if(verbose) message("Applying the Poisson mixture model...")

        # ---- block processing ----
        #block_size <- block_size
        blocks <- split(seq_len(n_cells),
                        ceiling(seq_len(n_cells) / block_size))
        #message(blocks[[1]])

        deleted_off_blocks <- rep(FALSE, if(use_offset_blocks) length(off_block_files) else 0L)

        for (b in seq_along(blocks)) {
          if(verbose) message("Processing block ", b, "/", length(blocks))
          block_results <- if(!is.null(clust)) {
            parLapply(cl = clust, blocks[[b]], worker_fun)
          } else {
            lapply(blocks[[b]], worker_fun)
          }
          if(verbose) message("Storing results for block ", b, "/", length(blocks))
          store_results(blocks[[b]], block_results)
          rm(block_results)

          if(verbose && !is.null(clust)){
            worker_mem <- denoist_get_worker_rss(clust)
            rss_mb <- worker_mem$rss_kb / 1024
            valid <- is.finite(rss_mb)
            if(any(valid)){
              ord <- order(rss_mb[valid], decreasing = TRUE)
              top_n <- min(3L, length(ord))
              top_idx <- which(valid)[ord[seq_len(top_n)]]
              top_workers <- paste(sprintf("pid=%d:%.1fMB", worker_mem$pids[top_idx], rss_mb[top_idx]), collapse = ", ")
              message(sprintf(
                "[debug] worker RSS after block %d/%d | min=%.1fMB median=%.1fMB max=%.1fMB | top: %s",
                b, length(blocks), min(rss_mb[valid]), stats::median(rss_mb[valid]), max(rss_mb[valid]), top_workers
              ))
            } else {
              message(sprintf("[debug] worker RSS after block %d/%d | unavailable", b, length(blocks)))
            }
          }

          if(use_offset_blocks){
            # Remove offset blocks whose column ranges are fully behind the current block.
            current_max_idx <- max(blocks[[b]])
            removable <- which(!deleted_off_blocks & off_block_ends <= current_max_idx)
            if(length(removable) > 0L){
              unlink(off_block_files[removable], force = TRUE)
              deleted_off_blocks[removable] <- TRUE
              if(verbose){
                message(sprintf("[debug] cleaned %d consumed offset block file(s); %d remaining",
                                length(removable), sum(!deleted_off_blocks)))
              }
            }

            # Drop worker-side cached block data before next iteration.
            if(!is.null(clust)){
              clusterEvalQ(clust, {
                assign(".denoist_offset_cache", NULL, envir = .denoist_cache_env)
                gc(FALSE)
                NULL
              })
            } else {
              assign(".denoist_offset_cache", NULL, envir = .denoist_cache_env)
            }
          }

          gc(FALSE)
        }

        if(!is.null(clust)){
          stopCluster(clust)
        }
    } else {
      if(verbose){
        message("Using in-memory processing...")
      }
      if(.Platform$OS.type == "unix" && cl > 1L){
        mc_cores <- max(1L, cl - 1L)
        if(verbose){
          message(sprintf("Using %d worker(s) with mclapply.", mc_cores))
          message("Applying the Poisson mixture model...")
        }
        results <- parallel::mclapply(
          seq_len(ncol(mat)),
          apply_poisson_mixture_single(get_count_col, get_offset_col, posterior_cutoff, n_inits, init_scheme, pmm_return_posterior, pmm_membership_output),
          mc.cores = mc_cores
        )
        store_results(seq_len(ncol(mat)), results)
        rm(results)
      } else {
        cluster_info <- denoist_init_cluster(cl - 1L, prefer_fork = FALSE, verbose = verbose)
        if(is.null(cluster_info)){
          warning("Cluster setup failed for all worker counts/types; falling back to sequential PMM execution.")
          if(verbose){
            message("Applying the Poisson mixture model without parallel processing...")
          }
          results <- pblapply(seq_len(ncol(mat)),
            apply_poisson_mixture_single(get_count_col, get_offset_col, posterior_cutoff, n_inits, init_scheme, pmm_return_posterior, pmm_membership_output),
                              cl = 1)
          store_results(seq_len(ncol(mat)), results)
          rm(results)
        } else {
          clust <- cluster_info$cluster
          if(verbose){
            message(sprintf("Using %d worker(s) with %s cluster.", cluster_info$workers, cluster_info$type))
            message("Applying the Poisson mixture model...")
          }
          results <- parLapply(cl = clust,
                  seq_len(ncol(mat)),
            apply_poisson_mixture_single(get_count_col, get_offset_col, posterior_cutoff, n_inits, init_scheme, pmm_return_posterior, pmm_membership_output))
          store_results(seq_len(ncol(mat)), results)
          rm(results)
          stopCluster(clust)
        }
      }
    }
  } else {
    # use pblapply normally without multiple processing
    if(verbose){
      message("Applying the Poisson mixture model without parallel processing...")
    }
    results <- pblapply(seq_len(ncol(mat)),
          apply_poisson_mixture_single(get_count_col, get_offset_col, posterior_cutoff, n_inits, init_scheme, pmm_return_posterior, pmm_membership_output),
              cl = cl)
    store_results(seq_len(ncol(mat)), results)
    rm(results)
  }

  # -------------------------------
  # Assemble results
  # -------------------------------
  if(verbose) message("Tidying up results...")

  if(use_sparse_memberships_out){
    memberships_out <- list(
      mode = "zero_idx",
      zero_idx_by_col = memberships_zero_idx_by_col,
      nrow = nrow(mat),
      ncol = ncol(mat),
      rownames = rownames(mat),
      colnames = colnames(mat)
    )
  } else {
    memberships_out <- memberships_matrix
  }

  if(use_sparse_memberships_out){
    adjusted_counts <- mat
    # Use block-wise parallelism to avoid dense membership reconstruction.
    if(cl > 1L && ncol(mat) > block_size && .Platform$OS.type == "unix"){
      if(verbose) message("Zeroing counts from sparse memberships using block-wise parallelism...")
      blocks <- split(seq_len(ncol(mat)), ceiling(seq_len(ncol(mat)) / block_size))
      block_worker <- function(cols){
        block <- mat[, cols, drop = FALSE]
        for(k in seq_along(cols)){
          idx <- cols[k]
          z <- memberships_zero_idx_by_col[[idx]]
          if(length(z) > 0L){
            block[z, k] <- 0L
          }
        }
        block
      }
      block_results <- pblapply(blocks, block_worker, cl = cl)
      adjusted_counts <- do.call(cbind, block_results)
      rm(block_results)
      if(verbose) message("Finished zeroing counts from sparse memberships.")
    } else {
      if(verbose) message("Zeroing counts from sparse memberships (sequential)...")
      for(idx in seq_len(ncol(mat))){ #bad idea to use a loop here but the alternative is to reconstruct a dense memberships matrix which we are trying to avoid by using the sparse zero_idx representation in the first place
        z <- memberships_zero_idx_by_col[[idx]]
        if(length(z) > 0L){
          adjusted_counts[z, idx] <- 0L
        }
      }
      if(verbose) message("Finished zeroing counts from sparse memberships.")
    }
  } else {
    adjusted_counts <- mat * memberships_matrix
  }
  colnames(adjusted_counts) <- colnames(mat)
  rownames(adjusted_counts) <- rownames(mat)

  # -------------------------------
  # Save results
  # -------------------------------
  if(!is.null(out_dir)){
    if(!dir.exists(out_dir)){
      dir.create(out_dir, recursive = TRUE)
    }

    saveRDS(list(
      memberships = memberships_out,
      adjusted_counts = adjusted_counts,
      params = params_out
    ), file = file.path(out_dir, "denoist_results.rds"))
  }

  return(list(
    memberships = memberships_out,
    adjusted_counts = adjusted_counts,
    params = params_out
  ))
}
