# Computes standard Kalman filter for linear, time-homogeneous Gaussian state space model, with initial state distribution known
# Notation is from D&K (see section 4.3), with A_t replacing T_t, and indices dropped for time-homogeneous version, like this:
#
# y_t = Z\alpha_t + \varepsilon_t, \varepsilon_t \sim N(0, H)
# \alpha_{t+1} = A\alpha_t + R\eta_t, \eta_t \sim N(0, Q)
# \alpha_1 \sim N(a1, P1)

kfilter <- function(yt, Z, A, R, H, Q, a1, P1) {
  m <- length(a1)
  n <- ifelse(is.matrix(yt), nrow(yt), length(yt))

  # Forecast distribution
  at <- matrix(nrow = n + 1, ncol = m)
  Pt <- array(dim = c(m, m, n + 1))

  # Filtered distribution
  att <- matrix(nrow = n, ncol = m)
  Ptt <- array(dim = c(m, m, n))

  # Initialization
  at[1, ] <- a1
  Pt[, , 1] <- P1

  for (t in 1:n) {
    vt <- yt[t] - Z %*% at[t, ]
    Ft <- Z %*% Pt[, , t] %*% t(Z) + H

    # Filtering step
    att[t, ] <- at[t, ] + Pt[, , t] %*% t(Z) %*% solve(Ft) %*% vt
    Ptt[, , t] <- Pt[, , t] - Pt[, , t] %*% t(Z) %*% solve(Ft) %*% Z %*% Pt[, , t]

    # Forecasting step
    at[t+1, ] <- A %*% att[t, ]
    Pt[, , t+1] <- A %*% Ptt[, , t] %*% t(A) + R %*% Q %*% t(R)
  }

  return(list(filtered.mean = att, filtered.variance = Ptt))
}
