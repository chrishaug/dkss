#' A linear Gaussian state-space model
#'
#' The model is described in D&K as follows:
#' \deqn{y_t = Z\alpha_t + \epsilon_t, \epsilon_t ~ N(0, H)}
#' \deqn{\alpha_{t+1} = A\alpha_t + R\eta_t, \eta_t ~ N(0, Q)}
#' \deqn{\alpha_1 ~ N(a_1, P_1)}
#'
#' @param yt The observed data, either as a vector or as a matrix with time as rows (n x p). Can also be a ts (mts) object.
#' @param Z (p x m) matrix, or a vector
#' @param A (m x m) matrix, or a scalar (called T in D&K, changed because of the TRUE boolean)
#' @param R (m x r) matrix, or a vector
#' @param H (p x p) matrix, or a scalar
#' @param Q (r x r) matrix, or a scalar
#' @param a1 (m x 1) matrix, or a vector
#' @param P1 (m x m) matrix, or a scalar
#'
#' @export
lgss <- function(yt, Z, A, R, H, Q, a1, P1) {
  # Need to convert any vectors to matrices to treat them uniformly
  yt <- as.matrix(yt)
  a1 <- as.matrix(a1)
  Q <- as.matrix(Q)

  # A few useful sizes
  m <- nrow(a1)
  n <- nrow(yt)
  p <- ncol(yt)
  r <- ncol(R)

  # For the ones that aren't matrices, this is our best guess at their size
  if (!is.matrix(Z)) dim(Z) <- c(p, m)
  if (!is.matrix(Z)) dim(A) <- c(m, m)
  if (!is.matrix(Z)) dim(R) <- c(m, r)
  if (!is.matrix(Z)) dim(H) <- c(p, p)
  if (!is.matrix(Z)) dim(P1) <- c(m, m)

  # Check if all the dimensions are coherent
  stopifnot(nrow(Z) == p, ncol(H) == p, nrow(H) == p)
  stopifnot(ncol(Z) == m, nrow(A) == m, ncol(A) == m, nrow(R) == m, nrow(P1) == m, ncol(P1) == m)
  stopifnot(nrow(Q) == r, ncol(Q) == r)

  result <- list(yt=yt, Z=Z, A=A, R=R, H=H, Q=Q, a1=a1, P1=P1)
  class(result) <- "lgss"

  return(result)
}

is.lgss <- function(mod) {
  return(inherits(mod, "lgss"))
}
