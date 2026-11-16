#' Fit an ordinary least squares regression model
#'
#' Fits a linear model using QR decomposition with modified Gram-Schmidt.
#' Rows with missing values in the model variables are omitted.
#' The design matrix must have full column rank.
#'
#' @param formula A model formula, such as `mpg ~ wt`.
#' @param data A data frame containing the variables in `formula`.
#' @return A list with three elements:
#'   - `coefficients`: The estimated regression coefficients.
#'   - `Q`: The matrix with orthonormal columns.
#'   - `R`: The upper triangular matrix.
#' @examples
#' # 1. One predictor: model fuel economy using vehicle weight.
#' fit_one <- matrix_lm(mpg ~ wt, data = mtcars)
#' fit_one$coefficients
#'
#' # 2. Multiple predictors: use weight and horsepower.
#' fit_multi <- matrix_lm(mpg ~ wt + hp, data = mtcars)
#' fit_multi$coefficients
#' all.equal(fit_multi$coefficients, coef(lm(mpg ~ wt + hp, data = mtcars)))
#'
#' # 3. Full column rank: reconstruct the design matrix from Q and R.
#' dat <- data.frame(y = c(2, 4, 5, 8), x = 0:3)
#' X <- model.matrix(y ~ x, data = dat)
#' qr(X)$rank == ncol(X)
#' fit_full <- matrix_lm(y ~ x, data = dat)
#' all.equal(unname(X), fit_full$Q %*% fit_full$R)
#'
#' # 4. Rank deficiency: x2 is exactly twice x.
#' dat$x2 <- 2 * dat$x
#' X_dep <- model.matrix(y ~ x + x2, data = dat)
#' qr(X_dep)$rank < ncol(X_dep)
#' tryCatch(
#'   matrix_lm(y ~ x + x2, data = dat),
#'   error = function(e) paste("Expected error:", conditionMessage(e))
#' )
#'
#' # 5. Missing predictor value: compare with explicitly dropping that row.
#' cars_na <- mtcars
#' cars_na$wt[1] <- NA_real_
#' fit_na <- matrix_lm(mpg ~ wt, data = cars_na)
#' fit_complete <- matrix_lm(mpg ~ wt, data = cars_na[-1, ])
#' all.equal(fit_na$coefficients, fit_complete$coefficients)
#'
#' # 6. Intercept-only model: estimate the mean fuel economy.
#' fit_intercept <- matrix_lm(mpg ~ 1, data = mtcars)
#' fit_intercept$coefficients
#' all.equal(unname(fit_intercept$coefficients), mean(mtcars$mpg))
#' @export
matrix_lm <- function(formula, data){
  #prepare and organize data
  dat <- model.frame(formula, data = data, na.action = na.omit)
  y <- model.response(dat) #outcome
  X <- model.matrix(terms(dat), data = dat) 
  
  #collect observations and predictors
  n <- nrow(X)
  p <- ncol(X)
  
  # check for full column rank
  d <- svd(X, nu = 0, nv = 0)$d
  cutoff <- 1e-10 * max(d)
  num_kept <- sum(d > cutoff)
  
  if (num_kept < p) {
    stop("Design matrix has numerical rank ", num_kept,
         " but has ", p, " columns.")
  }
  
  Q <- matrix(0, nrow = n, ncol = p)
  R <- matrix(0, nrow = p, ncol = p)
  beta <- numeric(p)
  
  for (j in seq_len(p)) {
    v <- X[, j]
    
    for (k in seq_len(j - 1)) {
      q <- Q[, k]
      R[k, j] <- sum(q * v)
      v <- v - R[k, j] * q
    }
    
    R[j, j] <- sqrt(sum(v^2))
    Q[, j] <- v / R[j, j]
  }
  
  qty_matrix <- t(Q) %*% y
  qty <- drop(qty_matrix)
  for (j in rev(seq_len(p))) {
    rhs <- qty[j]
    
    if (j < p) {
      right <- seq.int(j + 1L, p)
      rhs <- rhs - sum(R[j, right] * beta[right])
    }
    
    beta[j] <- rhs / R[j, j]
  }
  
  names(beta) <- colnames(X)
  list(coefficients = beta, Q = Q, R = R)
}

