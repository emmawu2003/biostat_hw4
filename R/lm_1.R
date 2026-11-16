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

