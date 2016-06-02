# Computes standard Kalman filter for linear, time-homogeneous Gaussian state space model, with initial state distribution known
# Notation is from D&K (see section 4.3), with A_t replacing T_t, and indices dropped for time-homogeneous version, like this:
#
# y_t = Z\alpha_t + \varepsilon_t, \varepsilon_t \sim N(0, H)
# \alpha_{t+1} = A\alpha_t + R\eta_t, \eta_t \sim N(0, Q)
# \alpha_1 \sim N(a1, P1)

kfilter <- function(mod) {
  # Check that mod is indeed of type lgss
  if (!is.lgss(mod)) stop("This function requires an object of type 'lgss'")

  # A few useful sizes
  m <- nrow(mod$a1)
  n <- nrow(mod$yt)
  p <- ncol(mod$yt)

  # Forecast distribution
  at <- matrix(nrow = m, ncol = n + 1)
  Pt <- array(dim = c(m, m, n + 1))
  vt <- matrix(nrow = p, ncol = n)
  Ft <- array(dim = c(p, p, n))

  # Filtered distribution
  att <- matrix(nrow = m, ncol = n)
  Ptt <- array(dim = c(m, m, n))

  # Initialization
  at[, 1] <- mod$a1
  Pt[, , 1] <- mod$P1

  for (t in 1:n) {
    vt[, t] <- mod$yt[t, ] - mod$Z %*% at[, t]
    Ft[, , t] <- mod$Z %*% Pt[, , t] %*% t(mod$Z) + mod$H

    # Filtering step
    att[, t] <- at[, t] + Pt[, , t] %*% t(mod$Z) %*% solve(Ft[, , t]) %*% vt[, t]
    Ptt[, , t] <- Pt[, , t] - Pt[, , t] %*% t(mod$Z) %*% solve(Ft[, , t]) %*% mod$Z %*% Pt[, , t]

    # Forecasting step
    at[, t+1] <- mod$A %*% att[, t]
    Pt[, , t+1] <- mod$A %*% Ptt[, , t] %*% t(mod$A) + mod$R %*% mod$Q %*% t(mod$R)
  }

  filtered <- lgss.filtered(att, Ptt, at[, 1:n], Pt[, , 1:n], vt, Ft)

  return(filtered)
}
