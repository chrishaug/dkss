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

is.lgss.smoothed <- function(x) {
  return(inherits(x, "lgss.smoothed"))
}
