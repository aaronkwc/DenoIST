test_mat <- readRDS(test_path("testdata", "test_mat.rds"))
test_coords <- readRDS(test_path("testdata", "test_coords.rds"))
test_tx <- readRDS(test_path("testdata", "test_tx.rds"))
test_spe <- readRDS(test_path("testdata", "test_spe.rds"))

test_that("DenoIST works with SpatialExperiment input", {
   # Test the function with the provided data
  res <- denoist(mat = test_spe,
              coords = NULL,
              tx = test_tx,
              distance = 50, nbins = 200, cl = 1,
              out_dir = "denoist_results", verbose = TRUE)
  expect_length(res, 3)
})


test_that("DenoIST works with Matrix input", {
  # Test the function with the provided data
  res <- denoist(mat = test_mat,
               coords = test_coords,
               tx = test_tx,
               distance = 50, nbins = 200, cl = 1, verbose = TRUE)
  expect_length(res, 3)
})

test_that("DenoIST can be sped up with less inits", {
  # Test the function with the provided data
  res <- denoist(mat = test_mat,
                 coords = test_coords,
                 tx = test_tx,
                 distance = 50, nbins = 200, cl = 1, n_inits = 5, verbose = TRUE)
  expect_length(res, 3)
})

test_that("DenoIST takes vector as init input", {
  i <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  # Test the function with the provided data
  res <- denoist(mat = test_mat,
                 coords = test_coords,
                 tx = test_tx,
                 distance = 50, nbins = 200, cl = 1, n_inits = i, verbose = TRUE)
  expect_length(res, 3)
})

test_that("DenoIST supports sparse PMM membership buffering", {
  res <- denoist(mat = test_mat,
                 coords = test_coords,
                 tx = test_tx,
                 distance = 50,
                 nbins = 200,
                 cl = 1,
                 sparse_membership = TRUE,
                 keep_posterior = FALSE,
                 verbose = TRUE)
  expect_length(res, 3)
  expect_type(res$memberships, "list")
  expect_equal(res$memberships$mode, "zero_idx")
})

test_that("DenoIST can return sparse memberships", {
  res_sparse <- denoist(mat = test_mat,
                        coords = test_coords,
                        tx = test_tx,
                        distance = 50,
                        nbins = 200,
                        cl = 1,
                        sparse_membership = TRUE,
                        verbose = TRUE)

  res_dense <- denoist(mat = test_mat,
                       coords = test_coords,
                       tx = test_tx,
                       distance = 50,
                       nbins = 200,
                       cl = 1,
                       sparse_membership = FALSE,
                       verbose = TRUE)

  expect_type(res_sparse$memberships, "list")
  expect_equal(res_sparse$memberships$mode, "zero_idx")
  expect_length(res_sparse$memberships$zero_idx_by_col, ncol(test_mat))

  reconstructed <- matrix(1,
                          nrow = res_sparse$memberships$nrow,
                          ncol = res_sparse$memberships$ncol,
                          dimnames = list(res_sparse$memberships$rownames,
                                          res_sparse$memberships$colnames))
  for(idx in seq_len(res_sparse$memberships$ncol)){
    z <- res_sparse$memberships$zero_idx_by_col[[idx]]
    if(length(z) > 0L){
      reconstructed[z, idx] <- 0L
    }
  }

  expect_equal(reconstructed, res_dense$memberships)
})
