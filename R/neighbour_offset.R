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
#' @return A matrix of local background counts for each gene in each cell.
#' @importFrom pbapply pblapply
#' @importFrom hexbin hexbin
#' @importFrom dplyr group_by summarise
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom sparseMatrixStats rowSums2
#' @import dplyr
#' @import flexmix
#' @export
local_offset_distance_with_background <- function(mat,
                                                  tx,
                                                  coords,
                                                  tx_x = "x",
                                                  tx_y = "y",
                                                  feature_label = "gene",
                                                  distance = 50,
                                                  nbins = 200,
                                                  cl = 1) {
  #print(mat[1:5, 1:5])
  #print(coords[1:5, 1:2])
  message('Calculating global background...')
  # filter by qv20
  tx <- tx[tx[['qv']] >= 20,]
  #print(nrow(tx))
  #print(head(tx))
  # Create hexagonal bins
  hex_bins <- hexbin(tx[[tx_x]], tx[[tx_y]], xbins = nbins, IDs = T) # Adjust xbins for resolution

  x_range <- diff(range(tx[,tx_x]))
  hex_radius <- x_range / hex_bins@xbins / sqrt(3)

  # Calculate the area of each hexbin
  hex_area <- (3 * sqrt(3) / 2) * hex_radius^2

  # Assign each transcript to a hexbin using the `hexbin` object
  tx$hexbin_id <- hex_bins@cID  # Use the `cID` slot to get the cell IDs for each point
  tx$feature_name <- tx[,feature_label]

  # Group by hexbin and gene to count occurrences
  gene_bin_counts <- tx %>%
    group_by(hexbin_id, feature_name) %>%
    summarise(count = dplyr::n(), .groups = "drop")

  # Create a matrix of gene by bin
  gene_bin_matrix <- gene_bin_counts %>%
    pivot_wider(names_from = hexbin_id, values_from = count, values_fill = 0) %>%
    column_to_rownames(var = "feature_name")

  #print(rownames(mat)[1:5])
  gene_bin_matrix <- gene_bin_matrix[rownames(mat),]
  #print(gene_bin_matrix[1:5, 1:5])
  bin_total <-colSums(gene_bin_matrix)
  # Extract empty bins inferred from GMM
  # Fit a Gaussian Mixture Model to colSums(gene_bin_matrix)
  message("Running GMM...")
  #print(bin_total[1:5])
  #gmm <- Mclust(bin_total, G = 2)
  mo1 <- FLXMRglm(family = "gaussian")
  mo2 <- FLXMRglm(family = "gaussian")
  flexfit <- flexmix(x ~ 1, data = data.frame(x=bin_total), k = 2, model = list(mo1, mo2))
  # Get the parameters of the GMM
  c1 <- parameters(flexfit, component=1)[[1]]
  c2 <- parameters(flexfit, component=2)[[1]]
  # Print the summary of the GMM
  # print(summary(gmm))

  # Identify the component with the smaller mean
  #gmm_means <- gmm$parameters$mean
  gmm_means <- c(c1[1], c2[1])
  smaller_mean_component <- which.min(gmm_means)

  empty_bin_matrix <- gene_bin_matrix[,clusters(flexfit) == smaller_mean_component]
  empty_bin_matrix <- empty_bin_matrix[,colSums(empty_bin_matrix) > 0]

  per_unit_sum <- rowSums(empty_bin_matrix)/(ncol(empty_bin_matrix) * hex_area)
  scaled_sum <- per_unit_sum * distance^2 * pi

  bg_offset <- ifelse(scaled_sum == 0, 1, ceiling(scaled_sum))


  # for each obs, get neighbours within distance
  # and then get the total count
  get_neighbors_within_distance <- function(coords, distance) {
    coords_mat <- as.matrix(coords)
    neighbors <- vector("list", nrow(coords))
    neighbors <- pblapply(seq_len(nrow(coords)), function(i) {
      dists <- sqrt(rowSums2((coords_mat - coords_mat[i, ])^2))
      which(dists <= distance)
    }, cl = cl)
    return(neighbors)
  }
  message("Finding neighbours...")
  neighbors <- get_neighbors_within_distance(coords[, c("x", "y")], distance)
  #print(neighbors[[1]])

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

  message("Calculating local offset...")
  res <- pblapply(seq_len(ncol(mat)), get_local_offset, neighbors, mat, cl = cl)
  res_mat <- do.call(cbind, res)
  colnames(res_mat) <- colnames(mat)

  # add bg_offset to every column of res_mat
  res_mat <- sweep(res_mat, 1, bg_offset, "+")

  return(res_mat)
}
