#' @importFrom pbapply pblapply
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment assay
#' @export
#' @title DenoIST
#' @description
#' DenoIST (Denoising Image-based Spatial Transcriptomics) is a method for
#' identifying and removing contamination artefacts in
#' image-based single-cell transcriptomics (IST) data.
#' It uses a transposed Poisson mixture model to identify contamination.
#' @param mat A matrix of counts (genes x cells), or a SpatialExperiment object.
#' @param coords A data frame of coordinates (cells x 2).
#' @param tx A data frame of transcript coordinates (transcripts x 2).
#' @param distance The maximum distance to consider for local background estimation.
#' @param nbins The number of bins to use for hexagonal binning.
#' @param cl The number of cores to use for parallel processing.
#' @param out_dir The output directory to save the results.
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
#' mat <- matrix(rpois(1000, lambda = 5), nrow = 10, ncol = 100)
#' coords <- data.frame(x = runif(100), y = runif(100))
#' tx <- data.frame(x = runif(100), y = runif(100), qv = runif(100, 0, 1))
#' # Run DenoIST
#' result <- denoist(mat, coords, tx, distance = 50, nbins = 200, cl = 8, out_dir = "denoist_results")
#' # Check results
#' print(result$memberships)
#' print(result$adjusted_counts)
#' print(result$params)
denoist <- function(mat, tx, coords = NULL, distance = 50, nbins = 200, cl = 1, out_dir = NULL){
  # TODO:check input type
  if(is(mat, "SpatialExperiment")){
    coords <- spatialCoords(mat)
    mat <- assay(mat)
  }else if(is.null(coords)){
    stop("coords must be provided")
  }

  #print(coords[1:5,1:2])
  #print(colnames(coords))
  #print(mat[1:5, 1:5])
  #print(tx[1:5,])
  # calculate neighbour_offset
  message("Calculating neighbour offset...")
  off_mat <- local_offset_distance_with_background(mat = mat,
                                                   coords = coords,
                                                   tx = tx,
                                                   distance = distance,
                                                   nbins = nbins,
                                                   cl = cl)

  # Apply the Poisson mixture model
  message("Applying the Poisson mixture model...")
  results <- pblapply(1:ncol(mat),
                      apply_poisson_mixture_single,
                      mat,
                      off_mat,
                      cl = cl)

  # return neighbour_offset, adjusted_counts, posterior, params
  message("Tidying up results...")
  memberships_matrix <- do.call(cbind, lapply(results, function(res) res["memberships"]))
  memberships_matrix <- do.call(cbind, memberships_matrix)
  colnames(memberships_matrix) <- colnames(mat)
  rownames(memberships_matrix) <- rownames(mat)

  adjusted_counts <- mat * memberships_matrix
  colnames(adjusted_counts) <- colnames(mat)
  rownames(adjusted_counts) <- rownames(mat)

  # save the results
  if(!is.null(out_dir)){
    if(!dir.exists(out_dir)){
      dir.create(out_dir, recursive = TRUE)
    }
    saveRDS(list(memberships = memberships_matrix,
                 adjusted_counts = adjusted_counts,
                 params = results),
            file = paste0(out_dir, "/spade_results.rds"))
  }

  return(list(memberships = memberships_matrix,
              adjusted_counts = adjusted_counts,
              params = results))
}
