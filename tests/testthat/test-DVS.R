test_that("DVS returns a valid result", {
  set.seed(123)

  n <- 50
  p <- 20

  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("X", 1:p)

  y <- x[, 1] + x[, 2] + rnorm(n)

  result <- DVS(x, y, B = 10, Threshold = 5)

  expect_true(is.list(result))
  expect_true("selected" %in% names(result))
  expect_true(is.data.frame(result$selected))
  expect_true(all(c("Variable", "Selection_Frequency") %in%
                    names(result$selected)))
})
