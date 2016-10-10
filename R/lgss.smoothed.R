#' The result of running the Kalman smoother on an lgss.filtered object
#'
#' The results of the algorithm are defined in D&K as follows
#' \deqn{(\alpha_t|y_1, ..., y_n) ~ N(\hat{\alpha}_t, V_t)}
#'
#' We store the time dimension as the last one in all cases.
#'
#' @param filt The filtered model result, as defined by an object of class lgss.filtered
#' @param alphat Smoothed mean, an (m x n) matrix
#' @param Vt Smoothed variance, (m x m x n) array
#'
#' @export
lgss.smoothed <- function(filt, alphat, Vt) {
  # Type checking
  stopifnot(is.lgss.filtered(filt))
  stopifnot(is.matrix(alphat))
  stopifnot(is.array(Vt))
  stopifnot(length(dim(Vt))==3)

  # Check the sizes
  m <- nrow(alphat)
  n <- ncol(alphat)

  # Check if all the dimensions are coherent
  stopifnot(dim(Vt) == c(m, m, n))

  result <- list(filt=filt, alphat=alphat, Vt=Vt)
  class(result) <- "lgss.smoothed"

  return(result)
}

#' Checks whether an object is of class lgss.smoothed
#' @param x The object to be tested
#'
#' @export
is.lgss.smoothed <- function(x) {
  return(inherits(x, "lgss.smoothed"))
}

#' @export
plot.lgss.smoothed <- function(x, y, ...) {
  vert <- floor(sqrt(nrow(x$alphat)))
  horz <- ceiling(nrow(x$alphat)/vert)

  graphics::par(mfrow=c(vert,horz))
  for (i in 1:nrow(x$alphat)) {
    sds <- apply(x$Vt,3,function(P) sqrt(P[i,i]))
    dat <- cbind(x$alphat[i,], x$alphat[i,]-1.96*sds, x$alphat[i,]+1.96*sds)
    graphics::matplot(dat, col="blue", lty=c(1,2,2),type="l",ylab = "hat{alpha}",xlab="t")
  }
  graphics::par(mfrow=c(1,1))
}
