# Test solve_poisson_mixture
test_that("solve_poisson_mixture works with simple input", {
  x <- c(3, 5, 2, 8, 6)
  s <- c(1, 1, 1, 1, 1)
  result <- solve_poisson_mixture(x, s, verbose = FALSE)

  expect_type(result, "list")
  expect_true(all(result$memberships %in% c(0, 1)))
  expect_equal(length(result$memberships), length(x))
  expect_equal(length(result$posterior), length(x))
})

test_that("solve_poisson_mixture handles zero s values", {
  x <- c(3, 5, 2, 8, 6)
  s <- c(1, 0, 1, 0, 1)
  result <- solve_poisson_mixture(x, s, verbose = FALSE)

  expect_equal(length(result$memberships), length(x))
  expect_equal(result$memberships[which(s == 0)], rep(1, sum(s == 0)))
})

test_that("solve_poisson_mixture handles single init value as input", {
  x <- c(3, 5, 2, 8, 6)
  s <- c(1, 0, 1, 0, 1)
  result <- solve_poisson_mixture(x, s, n_inits = 5, verbose = FALSE)

  expect_equal(length(result$memberships), length(x))
  expect_equal(result$memberships[which(s == 0)], rep(1, sum(s == 0)))
})

test_that("solve_poisson_mixture handles vector init values as input", {
  x <- c(3, 5, 2, 8, 6)
  s <- c(1, 0, 1, 0, 1)
  i <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  result <- solve_poisson_mixture(x, s, n_inits = i, verbose = FALSE)

  expect_equal(length(result$memberships), length(x))
  expect_equal(result$memberships[which(s == 0)], rep(1, sum(s == 0)))
})

test_that("solve_poisson_mixture checks valid init values", {
  x <- c(3, 5, 2, 8, 6)
  s <- c(1, 0, 1, 0, 1)

  expect_error(solve_poisson_mixture(x, s, n_inits = c(-0.1, 0.2), verbose = FALSE))
  expect_error(solve_poisson_mixture(x, s, n_inits = c(0.1, 1.5), verbose = FALSE))
})

test_that("solve_poisson_mixture can skip posterior output", {
  x <- c(3, 5, 2, 8, 6)
  s <- c(1, 1, 1, 1, 1)
  result <- solve_poisson_mixture(x, s, return_posterior = FALSE, verbose = FALSE)

  expect_true(is.null(result$posterior))
  expect_equal(length(result$memberships), length(x))
})

test_that("solve_poisson_mixture can return sparse membership", {
  x <- c(3, 5, 2, 8, 6)
  s <- c(1, 1, 1, 1, 1)
  result <- solve_poisson_mixture(x, s,
                                  membership_output = "sparse",
                                  return_posterior = FALSE,
                                  verbose = FALSE)

  expect_true(is.null(result$memberships))
  expect_type(result$membership_zero_idx, "integer")
  expect_equal(result$n, length(x))
})
