#' Kalman filtering for linear Gaussian state-space model
#'
#' Applies classic Kalman filter algorithm to specified model, assuming that the parameters of the
#' initial distribution are known.
#'
#' @return  An object of class lgss.filtered which holds all the results of the algorithm.
#'
#' @param mod A linear Gaussian state space model as described by an object of class lgss.
#'
#' @export
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
