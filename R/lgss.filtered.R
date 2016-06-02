# An S3 class which represents the result of the Kalman filter for lgss object model

lgss.filtered <- function(att, Ptt, at, Pt, vt, Ft) {
  # Type checking
  stopifnot(is.matrix(att), is.matrix(at), is.matrix(vt))
  stopifnot(is.array(Ptt), is.array(Pt), is.array(Ft))
  stopifnot(length(dim(Ptt))==3, length(dim(Pt))==3, length(dim(Ft))==3)

  # Check the sizes
  m <- nrow(att)
  n <- ncol(att)
  p <- nrow(vt)

  # Check if all the dimensions are coherent
  stopifnot(dim(Ptt) == c(m, m, n), dim(Pt) == c(m, m, n), dim(Ft) == c(p, p, n))
  stopifnot(nrow(at) == m, ncol(at) == n, ncol(vt) == n)

  result <- list(att=att, Ptt=Ptt, at=at, Pt=Pt, vt=vt, Ft=Ft)
  class(result) <- "lgss.filtered"

  return(result)
}

is.lgss.filtered <- function(x) {
  return(inherits(x, "lgss.filtered"))
}

residuals.lgss.filtered <- function(filtered) {
  # Standardize the residuals (is this right?)
  filtered$vt/apply(filtered$Ft, function(x) sqrt(diag(x)))
}
