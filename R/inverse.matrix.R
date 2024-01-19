#' Pseudo-Inverse of Numeric Vectors and Matrices
#'
#' Compute the Moore-Penrose's pseudo inverse of numeric vectors and
#' matrices (possibly non-square).
#'
#' @param X numeric vector or (non-square) matrix.
#'
#' @param eps numeric, small positive real in the open \code{(0, 0.1)}; Or
#' \code{NULL} (see section Details). \code{eps} is a positive tolerance value
#' below which any singular value of a matrix is considered zero.
#'
#' @param ... additional argument passed to or from other methods. Currently not
#' used.
#'
#' @usage
#' inverse.matrix (X, eps = NULL, ...)
#'
#' inverse.vector (X, ...)
#'
#' @details
#' When argument \code{eps} is \code{NULL} (the default), it is set to
#' \code{.Machine$double.eps * max(NROW(X), NCOL(X)) * max(d)} where \code{d} is
#' the vector of singular values of \code{X}.
#'
#' \code{inverse.matrix} internally calls \code{inverse.vector} when \code{X} is
#' a vector (including when it is a matrix with only one row or one column).
#'
#' \code{inverse.vector} always turns its first argument into a vector.
#'
#' @return a numeric \code{'matrix'} if \code{object} is of class
#' \code{matrix}, and a numeric \code{'vector'} if \code{object} is of class
#' \code{vector}.
#'
#' The returned object also has attributes:
#'  * \code{'m'} which is the minimum of the numbers of rows (if any) and column (if any)
#'  of \code{object} (\code{m = 1} for any vector);
#'  * \code{'rank'} which is the matrix rank of \code{object} or its returned inverse
#'  (\code{rank = 1} for any vector).
#'
#' @seealso \link[bigpcor]{inverse} for the generic function.
#'
#' @aliases inverse.vector
#' @export inverse.matrix
#' @export inverse.vector
#'
#' @exportS3Method inverse matrix
#' @exportS3Method inverse vector
#' @import methods
#'
#' @note
#' The computations uses singular value decomposition as described in Section
#' \code{3.6.3} in \insertCite{petersen2012matrix;textual}{bigpcor} (see also
#' \insertCite{barnett1990matrices;textual}{bigpcor}). For rectangular matrices
#' (with two or more rows and columns), singular value decomposition is
#' performed using \link{svd}.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' #* Inverse of the vector [1, 2, 3]
#' inverse (1:3)
#'
#' #* Inverse of the row vector [1, 2, 3]
#' inverse (t(1:3))
#'
#' #* Inverse of the column vector [1, 2, 3]
#' inverse (t(t(1:3)))
#'
#' #* Inverse of a rectangular matrix and its transpose
#' A <- matrix(c(5, 2, 3, 2, 3, 5), nrow = 3, ncol = 2)
#' A
#' inverse (A)
#'
#' inverse (t(A))
#'
#' #* Because A is tall and of full column rank, inverse (A) times A is identity
#' inverse (A) %*% A
#'
#' #* Define B as:
#' B <- cbind(A, A[,1])
#' B
#'
#' #* Because B is not of full rank, inverse B times B is not identity
#' inverse (B) %*% B

inverse.matrix <- function (X, eps = NULL, ...) {
  stopifnot(is.numeric(X))

  # If an object of class vector
  if (is.vector(X))
    return(inverse.vector (X))

  # Otherwise:
  # Dimensions of X
  size <- dim(X)

  # Return error if X is a 3 or more dimensional array
  if (length(size) > 2) {
    stop("'X' must be a vector or a matrix object: 3 or more dimensional objects not handled.")
  }

  size1 <- size == 1

  # If X is a scalar
  if (all(size1)) {
    return(structure(1/X,
                     m = 1,
                     rank = 1))
  }

  # If X is a row or column vector
  if (any(size1)) {
    # Compute Moore-Penrose inverse of the vector
    out <- inverse.vector (X)

    if (size1[1]) {
      # If X is a row vector, the inverse is a column vector
      return(structure(as.matrix(out),
                       m = 1,
                       rank = 1))
    }

    # If X is a column vector, the inverse is a row vector
    return(structure(t(out),
                     m = 1,
                     rank = 1))

  }

  # Check 'eps'
  stopifnot(is.numeric(eps) | is.null(eps))

  # Compute Moore-Penrose inverse of the matrix
  return(mp.inverse (X = X, eps = eps))
}

# This function REMOVES NAs from 'X' before inverting the vector
inverse.vector <- function (X, ...) {
  stopifnot(is.numeric(X))
  X <- c(X)

  # If X is the null vector or all NA or NaNs, anything can be its inverse
  if (all(is.na(X)) | all(X == 0))
    return(structure(rep(NaN, length.out = length(X)) ,
                     m = 1,
                     rank = 1))

  # Compute Moore-Penrose inverse, remove NA and NaNs first
  inv <- X / sum(X * X, na.rm = TRUE)

  return(structure(inv,
                   m = 1,
                   rank = 1))
}

# Moore-Penrose Generalized Inverse
mp.inverse <- function (X, eps = NULL) {
  stopifnot(is.numeric(X))
  X <- as.matrix(X)

  # Perform singular value decomposition
  s <- svd(X)
  d <- s$d

  # Number of singular values (min(dim(X)))
  m <- length(d)

  # Check 'eps'
  if (missing(eps) | is.null(eps))
    eps <- .Machine$double.eps * max(NROW(X), NCOL(X)) * max(d)
  else {
    stopifnot(is.numeric(eps) | is.null(eps))
    stopifnot(eps > 0, eps < 0.1)
  }

  # remove singular values ~ zero
  d <- d[d > eps]
  n <- length(d)

  # Inverse of positive singular values
  inv <- if (n > 1) diag(1/d) else 1/d

  # Add rows, columns and rows of zeros if X has null singular values
  if (n != m) {
    inv <- cbind(inv, matrix(0, nrow=n, ncol=(m - n)))
    inv <- rbind(inv, matrix(0, nrow=(m-n), ncol=m))
  }

  # compute the Moore-Penrose inverse
  inv <- s$v %*% inv %*% t(s$u)

  # set very small values to zero
  inv[abs(inv) < eps] <- 0


  return(structure(inv, m = m, rank = n))
}

mpinv <- mp.inverse

cholinv <- function(X, ...) {
  chol2inv(chol(X, ...))
}
