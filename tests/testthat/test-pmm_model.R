# Test solve_poisson_mixture
test_that("solve_poisson_mixture works with simple input", {
  x <- c(3, 5, 2, 8, 6)
  s <- c(1, 1, 1, 1, 1)
  result <- solve_poisson_mixture(x, s, verbose = FALSE)

  expect_type(result, "list")
  expect_true(all(result$memberships %in% c(0, 1)))
  expect_equal(length(result$memberships), length(x))
  expect_equal(length(result$posterior), length(x))
  expect_true(result$lambda1 > 0)
  expect_true(result$lambda2 > 0)
  expect_true(result$pi >= 0 && result$pi <= 1)
})

test_that("solve_poisson_mixture handles zero s values", {
  x <- c(3, 5, 2, 8, 6)
  s <- c(1, 0, 1, 0, 1)
  result <- solve_poisson_mixture(x, s, verbose = FALSE)

  expect_equal(length(result$memberships), length(x))
  expect_equal(result$memberships[which(s == 0)], rep(1, sum(s == 0)))
})

# Test apply_poisson_mixture_single
test_that("apply_poisson_mixture_single handles errors gracefully", {
  c_matrix <- matrix(c(3, 5, 2, 8, 6), ncol = 1)
  s_matrix <- matrix(c(1, 0, 1, 0, 1), ncol = 1)
  result <- apply_poisson_mixture_single(1, c_matrix, s_matrix)

  expect_type(result, "list")
  expect_equal(length(result$memberships), nrow(c_matrix))
  expect_equal(length(result$posterior), nrow(c_matrix))
})
