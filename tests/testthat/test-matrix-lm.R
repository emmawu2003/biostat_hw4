test_that("coefficients agree with lm for multiple predictors", {
  fit <- matrix_lm(
    Fertility ~ Education + Agriculture,
    data = swiss
  )
  
  reference <- stats::lm(
    Fertility ~ Education + Agriculture,
    data = swiss
  )
  
  expect_equal(
    fit$coefficients,
    coef(reference),
    tolerance = 1e-8
  )
})