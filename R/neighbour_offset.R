#' @title Neighbourhood offset
#' @description This function calculates the local neighbourhood offset together with ambient background for each cell in a count matrix.
#' @param mat A count matrix with genes as rows and cells as columns.
#' @param tx A transcript dataframe with x, y coordinates and qv values.
#' @param coords A dataframe with x, y coordinates of each cell as separate columns.
#' @param tx_x Column name for the x coordinates in the transcripts dataframe.
#' @param tx_y Column name for the y coordinates in the transcripts dataframe.
#' @param feature_label Column name for the gene of each transcript in the transcripts dataframe.
#' @param distance The maximum distance to consider for local background estimation.
#' @param nbins The number of bins to use for hexagonal binning, used for calculating background transcript contamination.
#' @param cl The number of cores to use for parallel processing.
#' @param verbose Logical, if TRUE, print progress messages.
#' @return A matrix of local background counts for each gene in each cell.
#' @details
#' The function calculates the offset used for each cell based on their local
#' neighbourhoods.In most cases you do not need to use this as denoist already
#' runs this internally but it is good for debugging if needed.
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
#' off_mat <- local_offset_distance_with_background(mat, tx, coords,
#'                                                  distance = 1, nbins = 50,
#'                                                  cl = 1, verbose = TRUE)
#' # Check results
#' print(off_mat[1:5, 1:5])
#' @importFrom pbapply pblapply
#' @importFrom hexbin hexbin
#' @importFrom sparseMatrixStats rowSums2 colSums2
#' @importFrom stats xtabs
#' @importFrom flexmix FLXMRglm flexmix parameters clusters
#' @importFrom dbscan frNN
#' @export
local_offset_distance_with_background <- function(mat,
                                                  tx,
                                                  coords,
                                                  tx_x = "x",
                                                  tx_y = "y",
                                                  feature_label = "gene",
                                                  distance = 50,
                                                  nbins = 200,
                                                  cl = 1,
                                                  verbose = FALSE) {
  log_mem <- function(step_label) {
    if (verbose) {
      gc_stats <- gc()
      mb_cols <- which(grepl("Mb", colnames(gc_stats), fixed = TRUE))
      if (length(mb_cols) > 0) {
        used_mb <- sum(gc_stats[, mb_cols[1]])
        message(sprintf("[debug] %s | approx used memory: %.1f MB", step_label, used_mb))
      } else {
        message(sprintf("[debug] %s | gc() returned no MB column", step_label))
      }
    }
  }

  # store the gene names first in case a gene gets wiped out because of QV
  gene_names_tx <- unique(tx[[feature_label]])
  log_mem("Initialized inputs")

  # filter by qv20 if the column exists (ie Xenium)
  if('qv' %in% colnames(tx)){
    message('QV column found! Filtering qv for high quality transcripts...')
    tx <- tx[tx[['qv']] >= 20,]
  }else{
    message('QV column not found! Skipping filtering...')
  }
  if(verbose){
    message(sprintf("[debug] tx rows after qv filtering: %d", nrow(tx)))
  }
  log_mem("After qv filtering")

  if(verbose){
    message('Calculating global background...')
  }

  # Create hexagonal bins
  hex_bins <- hexbin(tx[[tx_x]], tx[[tx_y]], xbins = nbins, IDs = TRUE) # Adjust xbins for resolution
  log_mem("After hexbin creation")

  x_range <- diff(range(tx[[tx_x]]))
  hex_radius <- x_range / hex_bins@xbins / sqrt(3)

  # Calculate the area of each hexbin
  hex_area <- (3 * sqrt(3) / 2) * hex_radius^2

  # Build the contingency matrix from lightweight vectors to avoid growing tx.
  hexbin_id <- hex_bins@cID
  feature_name <- tx[[feature_label]]
  gene_bin_matrix <- xtabs(~ feature_name + hexbin_id, sparse = TRUE)
  if(verbose){
    message(sprintf("[debug] gene_bin_matrix dims: %d genes x %d bins", nrow(gene_bin_matrix), ncol(gene_bin_matrix)))
  }
  log_mem("After gene_bin_matrix construction")

  if(nrow(gene_bin_matrix) < length(gene_names_tx)){
    message("Warning: Some genes in the original transcript data were not included in the hexbin matrix. This may be due to low quality transcripts being filtered out. Consider adjusting the QV threshold or checking the input data.")
    # Allocate once and fill existing rows to avoid padding + rbind copies.
    missing_genes <- setdiff(gene_names_tx, rownames(gene_bin_matrix))
    if(length(missing_genes) > 0){
      full_gene_bin_matrix <- matrix(0,
                                     nrow = length(gene_names_tx),
                                     ncol = ncol(gene_bin_matrix),
                                     dimnames = list(gene_names_tx, colnames(gene_bin_matrix)))
      full_gene_bin_matrix[rownames(gene_bin_matrix), ] <- gene_bin_matrix
      gene_bin_matrix <- full_gene_bin_matrix
      rm(full_gene_bin_matrix)
    }
  }
  log_mem("After missing-gene handling")

  # Keep the map so bg_offset can be aligned back to all genes in mat.
  idx_map <- match(rownames(mat), rownames(gene_bin_matrix))
  matched_genes <- !is.na(idx_map)
  idx <- idx_map[matched_genes]

  gene_bin_matrix <- gene_bin_matrix[idx, , drop = FALSE]
  if(anyNA(gene_bin_matrix)){
    gene_bin_matrix[is.na(gene_bin_matrix)] <- 0
  }
  if(verbose){
    message(sprintf("[debug] gene_bin_matrix matched rows: %d", nrow(gene_bin_matrix)))
  }
  log_mem("After gene_bin_matrix alignment to mat")

  bin_total <- colSums2(gene_bin_matrix)
  log_mem("After bin_total computation")

  # Extract empty bins inferred from GMM
  # Fit a Gaussian Mixture Model to colSums(gene_bin_matrix)
  if(verbose){
    message("Running GMM...")
  }
  mo1 <- FLXMRglm(family = "gaussian")
  mo2 <- FLXMRglm(family = "gaussian")
  bg_offset <- tryCatch(
        { flexfit <- flexmix(x ~ 1, data = data.frame(x=bin_total), k = 2, model = list(mo1, mo2))
          # Get the parameters of the GMM
          c1 <- parameters(flexfit, component=1)[[1]]
          c2 <- parameters(flexfit, component=2)[[1]]
          # Print the summary of the GMM

          # Identify the component with the smaller mean
          gmm_means <- c(c1[1], c2[1])
          smaller_mean_component <- which.min(gmm_means)

          empty_bin_matrix <- gene_bin_matrix[,clusters(flexfit) == smaller_mean_component]
          empty_bin_matrix <- empty_bin_matrix[,colSums(empty_bin_matrix) > 0]

          per_unit_sum <- rowSums(empty_bin_matrix)/(ncol(empty_bin_matrix) * hex_area)
          scaled_sum <- per_unit_sum * distance^2 * pi

          bg_offset <- ifelse(scaled_sum == 0, 1, ceiling(scaled_sum))
          bg_offset
        }, error = function(e){
          message("flexmix failed during GMM fit: ", e$message)
          message("Setting global background contamination to 1...")
          bg_offset <- rep(1, nrow(gene_bin_matrix))
          bg_offset
        }
  )

  bg_offset <- as.numeric(bg_offset)
  bg_offset_aligned <- rep(1, nrow(mat))
  bg_offset_aligned[matched_genes] <- bg_offset
  bg_offset <- bg_offset_aligned
  if(verbose){
    message(sprintf("[debug] bg_offset length: %d", length(bg_offset)))
  }
  log_mem("After bg_offset alignment")

  # Release transcript-derived intermediates before neighbour/PMM stages.
  rm_list <- c(
    "tx", "gene_names_tx", "hex_bins", "hexbin_id", "feature_name",
    "gene_bin_matrix", "missing_genes", "idx_map", "matched_genes",
    "idx", "bin_total", "bg_offset_aligned"
  )
  rm(list = intersect(rm_list, ls()), inherits = FALSE)
  gc(FALSE)

  # for each obs, get neighbours within distance
  # and then get the total count
  get_neighbors_within_distance <- function(coords, distance) {
    coords_mat <- as.matrix(coords)
    neighbors <- vector("list", nrow(coords))
    neighbors <- pblapply(seq_len(nrow(coords)), function(i) {
      diff <- sweep(coords_mat, 2, coords_mat[i, ], FUN = "-")
      dists <- sqrt(rowSums2(diff^2))
      which(dists <= distance)
    }, cl = cl)

    #nn = frNN(coords_mat,eps=distance)

    #neighbors <- nn$id
    return(neighbors)
  }

  if(verbose){
    message("Finding neighbours...")
  }

  neighbors <- get_neighbors_within_distance(coords[, c(1,2)], distance)
  if(verbose){
    neighbor_counts <- lengths(neighbors)
    message(sprintf("[debug] neighbours summary: min=%d median=%d max=%d", min(neighbor_counts), as.integer(stats::median(neighbor_counts)), max(neighbor_counts)))
  }
  log_mem("After neighbor search")

  get_local_offset <- function(idx, neighbors, mat) {
    if (length(neighbors[[idx]]) == 0) {
      offset <- rep(0, nrow(mat)) + mat[, idx]
    } else {
      if (length(neighbors[[idx]]) == 1) {
        offset <- mat[, neighbors[[idx]]] + mat[, idx]
      } else {
        offset <- rowSums2(mat[, neighbors[[idx]]]) + mat[, idx]
      }
    }
    return(offset)
  }

  if(verbose){
    message("Calculating local offset...")
  }
  res <- pblapply(seq_len(ncol(mat)), get_local_offset, neighbors, mat, cl = cl)
  res_mat <- matrix(0,
                    nrow = nrow(mat),
                    ncol = ncol(mat),
                    dimnames = list(rownames(mat), colnames(mat)))
  for(idx in seq_along(res)){
    res_mat[, idx] <- res[[idx]]
  }
  rm(res)
  log_mem("After local offset matrix assembly")

  # add bg_offset to every column of res_mat
  res_mat <- sweep(res_mat, 1, bg_offset, "+")
  log_mem("After adding bg_offset")

  return(res_mat)
}

#' @title Neighbourhood offset fast
#' @description This function calculates the local neighbourhood offset together with ambient background for each cell in a count matrix.
#' This is a faster version of the local_offset_distance_with_background function that uses dbscan for neighbor finding.
#' Code credit to Sam Zimmerman with minor modifications.
#' @param mat A count matrix with genes as rows and cells as columns.
#' @param tx A transcript dataframe with x, y coordinates and qv values.
#' @param coords A dataframe with x, y coordinates of each cell as separate columns.
#' @param tx_x Column name for the x coordinates in the transcripts dataframe.
#' @param tx_y Column name for the y coordinates in the transcripts dataframe.
#' @param feature_label Column name for the gene of each transcript in the transcripts dataframe.
#' @param distance The maximum distance to consider for local background estimation.
#' @param nbins The number of bins to use for hexagonal binning, used for calculating background transcript contamination.
#' @param on_disk Logical, if TRUE write output in column blocks to disk and return block metadata instead of an in-memory dense matrix.
#' @param on_disk_dir Optional directory used when on_disk is TRUE. If NULL, a temporary directory is created.
#' @param block_size Number of columns to process per block in the fast matrix multiplication step.
#' @param cl The number of cores to use for parallel processing.
#' @param verbose Logical, if TRUE, print progress messages.
#' @return A matrix of local background counts for each gene in each cell.
#' @details
#' The function calculates the offset used for each cell based on their local
#' neighbourhoods.In most cases you do not need to use this as denoist already
#' runs this internally but it is good for debugging if needed. Usage is the same as the non-fast version, will replace the original function in denoist in the future if this one is stable and faster.
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
#' off_mat <- local_offset_distance_with_background_fast(mat, tx, coords,
#'                                                  distance = 1, nbins = 50,
#'                                                  cl = 1, verbose = TRUE)
#' # Check results
#' print(off_mat[1:5, 1:5])
#' @importFrom pbapply pblapply
#' @importFrom hexbin hexbin
#' @importFrom sparseMatrixStats rowSums2
#' @importFrom stats xtabs
#' @importFrom flexmix FLXMRglm flexmix parameters clusters
#' @importFrom dbscan frNN
#' @importFrom Matrix sparseMatrix
#' @importFrom methods as
#' @export
local_offset_distance_with_background_fast <- function(mat,
                                                        tx,
                                                        coords,
                                                        tx_x = "x",
                                                        tx_y = "y",
                                                        feature_label = "gene",
                                                        distance = 50,
                                                        nbins = 200,
                                                        cl = 1,
                                                        verbose = FALSE,
                                                        on_disk = FALSE,
                                                        on_disk_dir = NULL,
                                                        block_size = 5000L) {
  log_mem <- function(step_label) {
    if (verbose) {
      gc_stats <- gc()
      mb_cols <- which(grepl("Mb", colnames(gc_stats), fixed = TRUE))
      if (length(mb_cols) > 0) {
        used_mb <- sum(gc_stats[, mb_cols[1]])
        message(sprintf("[debug] %s | approx used memory: %.1f MB", step_label, used_mb))
      } else {
        message(sprintf("[debug] %s | gc() returned no MB column", step_label))
      }
    }
  }

  # store the gene names first in case a gene gets wiped out because of QV
  gene_names_tx <- unique(tx[[feature_label]])
  log_mem("Initialized inputs")

  # filter by qv20 if the column exists (ie Xenium)
  if('qv' %in% colnames(tx)){
    message('QV column found! Filtering qv for high quality transcripts...')
    tx <- tx[tx[['qv']] >= 20,]
  }else{
    message('QV column not found! Skipping filtering...')
  }
  if(verbose){
    message(sprintf("[debug] tx rows after qv filtering: %d", nrow(tx)))
  }
  log_mem("After qv filtering")

  if(verbose){
    message('Calculating global background...')
  }

  # Create hexagonal bins
  hex_bins <- hexbin(tx[[tx_x]], tx[[tx_y]], xbins = nbins, IDs = TRUE) # Adjust xbins for resolution
  log_mem("After hexbin creation")

  x_range <- diff(range(tx[[tx_x]]))
  hex_radius <- x_range / hex_bins@xbins / sqrt(3)

  # Calculate the area of each hexbin
  hex_area <- (3 * sqrt(3) / 2) * hex_radius^2

  # Build the contingency matrix from lightweight vectors to avoid growing tx.
  hexbin_id <- hex_bins@cID
  feature_name <- tx[[feature_label]]
  gene_bin_matrix <- xtabs(~ feature_name + hexbin_id)
  if(verbose){
    message(sprintf("[debug] gene_bin_matrix dims: %d genes x %d bins", nrow(gene_bin_matrix), ncol(gene_bin_matrix)))
  }
  log_mem("After gene_bin_matrix construction")

  if(nrow(gene_bin_matrix) < length(gene_names_tx)){
    message("Warning: Some genes in the original transcript data were not included in the hexbin matrix. This may be due to low quality transcripts being filtered out. Consider adjusting the QV threshold or checking the input data.")
    # Allocate once and fill existing rows to avoid padding + rbind copies.
    missing_genes <- setdiff(gene_names_tx, rownames(gene_bin_matrix))
    if(length(missing_genes) > 0){
      full_gene_bin_matrix <- matrix(0,
                                     nrow = length(gene_names_tx),
                                     ncol = ncol(gene_bin_matrix),
                                     dimnames = list(gene_names_tx, colnames(gene_bin_matrix)))
      full_gene_bin_matrix[rownames(gene_bin_matrix), ] <- gene_bin_matrix
      gene_bin_matrix <- full_gene_bin_matrix
      rm(full_gene_bin_matrix)
    }
  }
  log_mem("After missing-gene handling")

  # Keep the map so bg_offset can be aligned back to all genes in mat.
  idx_map <- match(rownames(mat), rownames(gene_bin_matrix))
  matched_genes <- !is.na(idx_map)
  idx <- idx_map[matched_genes]

  gene_bin_matrix <- gene_bin_matrix[idx, , drop = FALSE]
  if(anyNA(gene_bin_matrix)){
    gene_bin_matrix[is.na(gene_bin_matrix)] <- 0
  }
  if(verbose){
    message(sprintf("[debug] gene_bin_matrix matched rows: %d", nrow(gene_bin_matrix)))
  }
  log_mem("After gene_bin_matrix alignment to mat")

  bin_total <- colSums(gene_bin_matrix)
  log_mem("After bin_total computation")

  # Extract empty bins inferred from GMM
  # Fit a Gaussian Mixture Model to colSums(gene_bin_matrix)
  if(verbose){
    message("Running GMM...")
  }
  mo1 <- FLXMRglm(family = "gaussian")
  mo2 <- FLXMRglm(family = "gaussian")
  bg_offset <- tryCatch(
    { flexfit <- flexmix(x ~ 1, data = data.frame(x=bin_total), k = 2, model = list(mo1, mo2))
    # Get the parameters of the GMM
    c1 <- parameters(flexfit, component=1)[[1]]
    c2 <- parameters(flexfit, component=2)[[1]]
    # Print the summary of the GMM

    # Identify the component with the smaller mean
    gmm_means <- c(c1[1], c2[1])
    smaller_mean_component <- which.min(gmm_means)

    empty_bin_matrix <- gene_bin_matrix[,clusters(flexfit) == smaller_mean_component]
    empty_bin_matrix <- empty_bin_matrix[,colSums(empty_bin_matrix) > 0]

    per_unit_sum <- rowSums(empty_bin_matrix)/(ncol(empty_bin_matrix) * hex_area)
    scaled_sum <- per_unit_sum * distance^2 * pi

    bg_offset <- ifelse(scaled_sum == 0, 1, ceiling(scaled_sum))
    bg_offset
    }, error = function(e){
      message("flexmix failed during GMM fit: ", e$message)
      message("Setting global background contamination to 1...")
      bg_offset <- rep(1, nrow(gene_bin_matrix))
      bg_offset
    }
  )

  bg_offset <- as.numeric(bg_offset)
  bg_offset_aligned <- rep(1, nrow(mat))
  bg_offset_aligned[matched_genes] <- bg_offset
  bg_offset <- bg_offset_aligned
  if(verbose){
    message(sprintf("[debug] bg_offset length: %d", length(bg_offset)))
  }
  log_mem("After bg_offset alignment")

  # Release transcript-derived intermediates before neighbour/PMM stages.
  rm_list <- c(
    "tx", "gene_names_tx", "hex_bins", "hexbin_id", "feature_name",
    "gene_bin_matrix", "missing_genes", "idx_map", "matched_genes",
    "idx", "bin_total", "bg_offset_aligned"
  )
  rm(list = intersect(rm_list, ls()), inherits = FALSE)
  gc(FALSE)

  # for each obs, get neighbours within distance
  # and then get the total count
  get_neighbors_within_distance <- function(coords, distance) {
    coords_mat <- as.matrix(coords)
    nn = frNN(coords_mat,eps=distance)

    neighbors <- nn$id
    return(neighbors)
  }

  if(verbose){
    message("Finding neighbours...")
  }
  # note, self matches are not returned
  neighbors <- get_neighbors_within_distance(coords[, c(1,2)], distance)
  if(verbose){
    neighbor_counts <- lengths(neighbors)
    message(sprintf("[debug] neighbours summary: min=%d median=%d max=%d", min(neighbor_counts), as.integer(stats::median(neighbor_counts)), max(neighbor_counts)))
  }
  log_mem("After neighbor search")

  if(verbose){
    message("Calculating local offset...")
  }

  # alternatively try and combine all this into a single matrix algebra step
  n_cells <- ncol(mat)
  # 1. Create vectors of the row (i) and column (j) coordinates for our adjacency matrix
  # 'j' is the target cell index, repeated for as many neighbors as it has
  j_indices <- rep(seq_len(n_cells), lengths(neighbors))
  i_indices <- unlist(neighbors)

  A <- sparseMatrix(i = i_indices,
                    j = j_indices,
                    x = 1,
                    dims = c(n_cells, n_cells))
  A <- as(A > 0, "dMatrix")
  diag(A) <- 1
  log_mem("After adjacency matrix build")

  if(on_disk){
    block_size <- as.integer(block_size)
    if(is.na(block_size) || block_size < 1L){
      stop("block_size must be a positive integer.")
    }
    block_starts <- seq.int(1L, ncol(mat), by = block_size)
    if(is.null(on_disk_dir)){
      on_disk_dir <- tempdir()
    }
    dir.create(on_disk_dir, recursive = TRUE, showWarnings = FALSE)
    on_disk_run_dir <- tempfile("local_offset_run_", tmpdir = on_disk_dir)
    dir.create(on_disk_run_dir, recursive = TRUE, showWarnings = FALSE)

    block_files <- character(length(block_starts))
    block_col_ranges <- matrix(0,
                               nrow = length(block_starts),
                               ncol = 2,
                               dimnames = list(NULL, c("start", "end")))
    for (block_idx in seq_along(block_starts)) {
      block_start <- block_starts[block_idx]
      block_end <- min(block_start + block_size - 1L, ncol(mat))
      block_cols <- block_start:block_end
      block_col_ranges[block_idx, ] <- c(block_start, block_end)

      block_res <- as.matrix(mat %*% A[, block_cols, drop = FALSE])
      if(nrow(block_res) != nrow(mat) || ncol(block_res) != length(block_cols)){
        block_res <- matrix(block_res,
                            nrow = nrow(mat),
                            ncol = length(block_cols))
      }
      mat_block <- as.matrix(mat[, block_cols, drop = FALSE])
      if(nrow(mat_block) != nrow(mat) || ncol(mat_block) != length(block_cols)){
        mat_block <- matrix(mat_block,
                            nrow = nrow(mat),
                            ncol = length(block_cols))
      }
      block_res <- block_res + mat_block
      block_res <- sweep(block_res, 1, bg_offset, "+")

      if(nrow(block_res) != nrow(mat) || ncol(block_res) != length(block_cols)){
        stop(sprintf("Block dimension mismatch at block %d: expected %d x %d, got %d x %d",
                     block_idx, nrow(mat), length(block_cols), nrow(block_res), ncol(block_res)))
      }

      block_file <- file.path(on_disk_run_dir, sprintf("offset_block_%05d.rds", block_idx))
      saveRDS(list(cols = block_cols,
                   colnames = colnames(mat)[block_cols],
                   data = block_res),
              file = block_file)
      block_files[block_idx] <- block_file

      if(verbose && (block_idx == 1L || block_idx %% 10L == 0L || block_idx == length(block_starts))){
        message(sprintf("[debug] processed block %d/%d (%d:%d)",
                        block_idx, length(block_starts), block_start, block_end))
        log_mem("During blockwise accumulation")
      }
    }
    log_mem("After final offset accumulation")

    if(verbose){
      message(sprintf("[debug] on-disk mode wrote %d blocks to: %s", length(block_files), on_disk_run_dir))
    }
    return(list(mode = "on_disk_blocks",
                directory = on_disk_run_dir,
                files = block_files,
                block_col_ranges = block_col_ranges,
                nrow = nrow(mat),
                ncol = ncol(mat),
                cleanup_on_exit = FALSE,
                rownames = rownames(mat),
                colnames = colnames(mat)))
  }

  res_mat <- as.matrix(mat %*% A)
  log_mem("After matrix multiplication")
  res_mat <- res_mat + mat
  res_mat <- sweep(res_mat, 1, bg_offset, "+")
  log_mem("After final offset accumulation")

  colnames(res_mat) <- colnames(mat)
  rownames(res_mat) <- rownames(mat)
  return(res_mat)
}

denoist_compute_bg_offset <- function(mat,
                                      tx,
                                      tx_x,
                                      tx_y,
                                      feature_label,
                                      distance,
                                      nbins,
                                      verbose = FALSE) {
  # store the gene names first in case a gene gets wiped out because of QV
  gene_names_tx <- unique(tx[[feature_label]])

  # filter by qv20 if the column exists (ie Xenium)
  if("qv" %in% colnames(tx)){
    tx <- tx[tx[["qv"]] >= 20,]
  }

  # Create hexagonal bins
  hex_bins <- hexbin(tx[[tx_x]], tx[[tx_y]], xbins = nbins, IDs = TRUE)

  x_range <- diff(range(tx[[tx_x]]))
  hex_radius <- x_range / hex_bins@xbins / sqrt(3)
  hex_area <- (3 * sqrt(3) / 2) * hex_radius^2

  hexbin_id <- hex_bins@cID
  feature_name <- tx[[feature_label]]
  gene_bin_matrix <- xtabs(~ feature_name + hexbin_id)

  if(nrow(gene_bin_matrix) < length(gene_names_tx)){
    missing_genes <- setdiff(gene_names_tx, rownames(gene_bin_matrix))
    if(length(missing_genes) > 0){
      full_gene_bin_matrix <- matrix(0,
                                     nrow = length(gene_names_tx),
                                     ncol = ncol(gene_bin_matrix),
                                     dimnames = list(gene_names_tx, colnames(gene_bin_matrix)))
      full_gene_bin_matrix[rownames(gene_bin_matrix), ] <- gene_bin_matrix
      gene_bin_matrix <- full_gene_bin_matrix
      rm(full_gene_bin_matrix)
    }
  }

  idx_map <- match(rownames(mat), rownames(gene_bin_matrix))
  matched_genes <- !is.na(idx_map)
  idx <- idx_map[matched_genes]

  gene_bin_matrix <- gene_bin_matrix[idx, , drop = FALSE]
  if(anyNA(gene_bin_matrix)){
    gene_bin_matrix[is.na(gene_bin_matrix)] <- 0
  }

  bin_total <- colSums(gene_bin_matrix)

  if(verbose){
    message("Running GMM...")
  }
  mo1 <- FLXMRglm(family = "gaussian")
  mo2 <- FLXMRglm(family = "gaussian")
  bg_offset <- tryCatch(
    {
      flexfit <- flexmix(x ~ 1, data = data.frame(x=bin_total), k = 2, model = list(mo1, mo2))
      c1 <- parameters(flexfit, component=1)[[1]]
      c2 <- parameters(flexfit, component=2)[[1]]
      gmm_means <- c(c1[1], c2[1])
      smaller_mean_component <- which.min(gmm_means)

      empty_bin_matrix <- gene_bin_matrix[, clusters(flexfit) == smaller_mean_component]
      empty_bin_matrix <- empty_bin_matrix[, colSums(empty_bin_matrix) > 0]

      per_unit_sum <- rowSums(empty_bin_matrix) / (ncol(empty_bin_matrix) * hex_area)
      scaled_sum <- per_unit_sum * distance^2 * pi

      ifelse(scaled_sum == 0, 1, ceiling(scaled_sum))
    },
    error = function(e){
      message("flexmix failed during GMM fit: ", e$message)
      rep(1, nrow(gene_bin_matrix))
    }
  )

  bg_offset <- as.numeric(bg_offset)
  bg_offset_aligned <- rep(1, nrow(mat))
  bg_offset_aligned[matched_genes] <- bg_offset
  bg_offset_aligned
}

denoist_get_neighbors_fast <- function(coords, distance) {
  coords_mat <- as.matrix(coords)
  nn <- frNN(coords_mat, eps = distance)
  nn$id
}

