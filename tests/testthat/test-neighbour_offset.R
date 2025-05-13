test_that("Neighbour offsets can be calculated", {
  # Test the function with the provided data
  test_tx <- readRDS(test_path("testdata", "test_tx.rds"))
  test_mat <- readRDS(test_path("testdata", "test_mat.rds"))
  test_coords <- readRDS(test_path("testdata", "test_coords.rds"))

  off <- local_offset_distance_with_background(mat = test_mat,
                                              coords = test_coords,
                                              tx = test_tx,
                                              distance = 50, nbins = 200, cl = 1)
  expect_equal(nrow(off), nrow(test_mat))
  expect_equal(ncol(off), ncol(test_mat))
})
