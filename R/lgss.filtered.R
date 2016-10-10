#' The result of running the Kalman filter on an lgss object
#'
#' The results of the algorithm are defined in D&K as follows
#' \deqn{(\alpha_t|y_1, ..., y_{t-1}) ~ N(at, Pt)}
#' \deqn{(\alpha_t|y_1, ..., y_t) ~ N(att, Ptt)}
#' \deqn{v_t = y_t - Z_t a_t}
#' \deqn{(\v_t|y_1, ..., y_{t-1}) ~ N(0, Ft)}
#'
#' We store the time dimension as the last one in all cases.
#'
#' @param mod The model, as defined by an object of class lgss
#' @param att Filtered mean, an (m x n) matrix
#' @param Ptt Filtered variance, (m x m x n) array
#' @param at Forecast mean, an (m x n) matrix
#' @param Pt Forecast variance, (m x m x n) array
#' @param vt One-step ahead forecast errors, (p x n) matrix
#' @param Ft One-step ahead forecast variance, (p x p x n) array
#'
#' @export
lgss.filtered <- function(mod, att, Ptt, at, Pt, vt, Ft) {
  # Type checking
  stopifnot(is.lgss(mod))
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

  result <- list(mod=mod, att=att, Ptt=Ptt, at=at, Pt=Pt, vt=vt, Ft=Ft)
  class(result) <- "lgss.filtered"

  return(result)
}

#' Checks whether an object is of class lgss.filtered
#' @param x The object to be tested
#'
#' @export
is.lgss.filtered <- function(x) {
  return(inherits(x, "lgss.filtered"))
}

#' @export
residuals.lgss.filtered <- function(object, ...) {
  object$vt
}

#' @export
plot.lgss.filtered <- function(x, y, ...) {
  vert <- floor(sqrt(nrow(x$att)))
  horz <- ceiling(nrow(x$att)/vert)

  graphics::par(mfrow=c(vert,horz))
  for (i in 1:nrow(x$att)) {
    sds <- apply(x$Ptt,3,function(P) sqrt(P[i,i]))
    dat <- cbind(x$att[i,], x$att[i,]-1.96*sds, x$att[i,]+1.96*sds)
    graphics::matplot(dat, col="red", lty=c(1,2,2),type="l",ylab = "alpha_{t|t}",xlab="t")
  }
  graphics::par(mfrow=c(1,1))
}
