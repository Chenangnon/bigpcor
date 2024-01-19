
#' Rank-1 Update of Inverses of Inner Products
#'
#' Update a square matrix \eqn{A}, inverse of the inner product matrix
#' \code{t(X) %*% X}, after \eqn{X} has been updated by including an
#' additional column \eqn{v} in \eqn{X}, or by removing one column from \eqn{X}.
#'
#' @param object numeric matrix \eqn{A} of size \eqn{p \times p}, inverse of
#' the inner product matrix \eqn{t(X) \times X}. \code{object} can have an
#' attribute 'X', which is required if \code{add = TRUE} and argument \eqn{X} is missing.
#'
#' @param add logical, adding a new column into \eqn{X}?
#'
#' @param position a positive integer, indicating the position \eqn{j} of:
#' \itemize{
#'  \item{the vector \eqn{v} being added into \eqn{X} when \code{add = TRUE}}{, here \eqn{1 \leq j \leq p + 1} with default value \code{position}=\eqn{j}=\code{p+1};}
#'  \item{the vector \eqn{v} being dropped from \eqn{X} when \code{add = FALSE}}{, here \eqn{1 \leq j \leq p} with default value \code{position}=\eqn{j}=\code{p}.}
#' }
#'
#' @param v a length \eqn{n} numeric vector, required if \code{add = TRUE} (the
#' default) and not required nor used when \code{add = FALSE}.
#'
#' @param X \eqn{n \times p} numeric matrix. Note that if \eqn{X} is the centered
#' version of another data matrix \eqn{Y} divided by the square root of the
#' sample degrees of freedom \eqn{n-1}, then \eqn{A} is the inverse of the
#' sample estimate of the covariance matrix of \eqn{Y}.
#'
#' \eqn{X} can be missing as an argument, it must then be given as a named
#' attribute of \code{object} if \code{add = TRUE}. \eqn{X} is not required nor
#' used when \code{add = FALSE}.
#'
#' @param MP logical or \code{NULL}, is \code{object} a Moore-Penrose inverse?
#' The default (\code{NULL}) looks for a \code{rank} attribute of \code{object}
#' and is \code{MP = TRUE} if this rank in less than the size of \code{object},
#' and \code{MP = FALSE} otherwise. If \code{object} has no \code{rank}
#' attribute, the rank of \code{object} is determined using QR decomposition
#' (\link{qr}).
#'
#' @param ... additional argument passed to or from other methods. Currently not
#' used.
#'
#' @return an object of class \code{"inner.inv"}, i.e. a matrix with a named
#' attribute \code{"X"} (which can be \code{NULL}, i.e. missing).
#'
#' @seealso \link[bigpcor]{inverse} for computing pseudo-inverse of vector
#' and (non-square) matrices.
#'
#' @export update.inner.inv
#'
#' @exportS3Method update inner.inv
#' @import methods
#'
#' @note
#' The formula for a square non-singular matrix are taken from
#' \insertCite{khan2008updating;textual}{bigpcor} (see also Section \code{3.2.6}
#' in \insertCite{petersen2012matrix;textual}{bigpcor}).
#'
#' The equivalent formula for pseudo-inverse was obtained by combining Eq. (168)
#' (Section \code{3.2.7}) in \insertCite{petersen2012matrix;textual}{bigpcor}
#' with Eq. (2.1) in \insertCite{baksalary2007particular;textual}{bigpcor}.
#'
#' @references
#' \insertAllCited{}
#'
update.inner.inv <- function (object, add = TRUE,
                              position = NULL,
                              v, X, MP = NULL, ...) {
  stopifnot(is.matrix(object))
  p <- NCOL(object)
  stopifnot(NROW(object) == p)

  if (add) {
    stopifnot(!missing(v))
    stopifnot(is.numeric(object))
    if (missing(X)) {
      if (inherits(object, "inner.inv")) {
        X <- attr(object, "X")
        if(is.null(X)) {
          stop ("'object' must be a proper 'inner.inv' object with an 'X' attribute when 'X' is missing")
        }
      }
      else {
        X <- attr(object, "X")
        if(is.null(X)) {
          stop("'X' must be specifed if 'object' does not have an attribute 'X'")
        }
      }
    }
    stopifnot(all(NCOL(X) == dim(object)))
    stopifnot(length(v) == NROW(X))

    ### update of the inverse of A = t(X) %*% X after
    #   updating X by adding the vector v to X as X <- [X, v]
    inv <- add.last.inner.inv (A = object, X = X, v = v, MP = MP[1])

    # Position the added column to the desired position if not yet
    if (!is.null(position)) {
      stopifnot(is.numeric(position), length(position) == 1)
      stopifnot(floor(position) == position, position > 0, position <= p + 1)

      # Use a permutation matrix to re-position the vector v as the last
      if (position < p + 1) {
        Pmat <- diag(p + 1)
        if (position == 1) {
          Pmat <- cbind(Pmat[, p + 1], Pmat[, -(p + 1)])

          attr(inv, "X") <- cbind(v, X)

        }
        else {
          Pmat <- cbind(Pmat[, 1:(position - 1)], Pmat[, p + 1], Pmat[, position:p])

          attr(inv, "X") <- cbind(X[, 1:(position - 1)], v, X[, position:p])

        }

        # Switch columns and rows in A = object; t(Pmat) stands for solve(Pmat)
        inv <- t(Pmat) %*% inv %*% Pmat
      }
    }
    else {
      attr(inv, "X") <- cbind(X, v)
    }

    return(inv)
  }

  # Position the column to be dropped as the last if not yet
  if (!is.null(position)) {
    stopifnot(is.numeric(position), length(position) == 1)
    stopifnot(floor(position) == position, position > 0, position <= p)

    # Use a permutation matrix to re-position the vector v as the last
    if (position < p) {
      Pmat <- diag(p)
      Pmat <- cbind(Pmat[,-position], Pmat[,position])

      # Switch columns and rows in A = object; t(Pmat) stands for solve(Pmat)
      object <- t(Pmat) %*% object %*% Pmat
    }
  }

  ### update of the inverse of A = t(X) %*% X after
  #   updating X by dropping the vector v from as: X = [X0, v] --> X = X0
  return(drop.last.inner.inv (A = object, MP = MP[1]))
}

# matupdate.inner.inv <- update.inner.inv

## Implement for MP = TRUE rankOut
add.last.inner.inv <- function (A, X, v, MP = FALSE) {

  p <- NCOL(X)
  stopifnot(c(NROW(A), NCOL(A)) == p)
  v <- as.vector(v)
  stopifnot(is.numeric(v), length(v) == NROW(X))
  stopifnot(is.numeric(MP) | is.logical(MP))
  MP <- as.logical(MP)
  MP <- MP[1]
  if (is.null(MP)) {
    MP <- attr(A, "rank")
    if (is.null(MP)) {
      MP <- qr(A)$rank < p
    }
    else {
      MP <- MP < p
    }
  }

  if (MP) { # If A is a Moore-Penrose inverse



    stop("not yet implemented")

    rankOut <- qr(cbind(X, v))$rank

    return(structure(A, rank = rankOut, class = c("matrix", "inner.inv")))
  }
  else { # Otherwise
    Xtv <- c(t(v) %*% X)
    AXtv <- c(A %*% Xtv)
    d <- sum(v*v) - sum(Xtv * AXtv)
    inv <- rbind(
      cbind(A * d + outer(AXtv, AXtv), - AXtv),
      c(- AXtv, 1) ) / d

    return(structure(inv, rank = p + (d != 0),
                     class = c("matrix", "inner.inv")))
  }
}

## Implement for MP = TRUE
drop.last.inner.inv <- function (A, MP = NULL) {
  p <- NCOL(A)
  stopifnot(NROW(A) == p)
  MP <- MP[1]
  if (is.null(MP)) {
    MP <- attr(A, "rank")
    if (is.null(MP)) {
      MP <- qr(A)$rank < p
    }
    else {
      MP <- MP < p
    }
  }

  if (MP) {
    stop("not yet implemented")


    return(structure(inv,
                     rank = p - 1,
                     class = c("matrix", "inner.inv")))
  }
  else {
    f <- A[1:(p - 1), p]
    inv <- A[1:(p - 1), 1:(p - 1)] - outer (f, f) / A[p,p]

    return(structure(inv,
                     rank = p - 1,   ### REVISIT THIS
                     class = c("matrix", "inner.inv")))
  }
}
