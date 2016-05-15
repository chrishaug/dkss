# Computes standard Kalman filter for linear, time-homogeneous Gaussian state space model, with initial state distribution known
# Notation is from D&K (see section 4.3), with A_t replacing T_t, and indices dropped for time-homogeneous version, like this:
#
# y_t = Z\alpha_t + \varepsilon_t, \varepsilon_t \sim N(0, H)
# \alpha_{t+1} = A\alpha_t + R\eta_t, \eta_t \sim N(0, Q)
# \alpha_1 \sim N(a1, P1)

kfilter <- function(yt, Z, A, R, H, Q, a1, P1) {
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

  # Forecast distribution
  at <- matrix(nrow = m, ncol = n + 1)
  Pt <- array(dim = c(m, m, n + 1))

  # Filtered distribution
  att <- matrix(nrow = m, ncol = n)
  Ptt <- array(dim = c(m, m, n))

  # Initialization
  at[, 1] <- a1
  Pt[, , 1] <- P1

  for (t in 1:n) {
    vt <- yt[t, ] - Z %*% at[, t]
    Ft <- Z %*% Pt[, , t] %*% t(Z) + H

    # Filtering step
    att[, t] <- at[, t] + Pt[, , t] %*% t(Z) %*% solve(Ft) %*% vt
    Ptt[, , t] <- Pt[, , t] - Pt[, , t] %*% t(Z) %*% solve(Ft) %*% Z %*% Pt[, , t]

    # Forecasting step
    at[, t+1] <- A %*% att[, t]
    Pt[, , t+1] <- A %*% Ptt[, , t] %*% t(A) + R %*% Q %*% t(R)
  }

  return(list(filtered.mean = att, filtered.variance = Ptt))
}
