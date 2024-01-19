#' Generalized Inverse of Numeric Objects
#'
#' \code{inverse} is a generic function whose default methods compute the
#' Moore-Penrose's generalized inverse of numeric vectors and matrices.
#'
#' @export inverse
#'
#' @param object a numeric matrix-like object of a class with a
#' \code{inverse} method. For the default method, \code{object} is a numeric
#' vector or a matrix.
#'
#' @param ... additional argument passed to or from other methods.
#'
#' Currently, the default method uses \code{eps}: a positive tolerance
#' below which any singular value of a matrix is considered zero.
#'
#' @exportMethod inverse
#'
#' @details
#' Additional \code{inverse} methods can be defined for specific classes.
#'
#' For the default method, when \code{object} is a numeric vector, or a matrix
#' of only one row or one column, any \code{NA} or \code{NaN} is ignored.
#' If \code{object} is a matrix with two rows or more and two columns or more,
#' any \code{NA} or \code{NaN} in \code{object} causes error.
#'
#' @return An object of the same class as \code{object}.
#'
#' For the default method
#' \link{inverse.matrix}, a numeric \code{'matrix'} if \code{object} is of class
#' \code{matrix}, and a numeric vector if \code{object} is of class \code{vector}.
#' The returned object also has attributes:
#'  * \code{m} the minimum of the number of rows (if any),  and column, if any,
#'  of \code{object};
#'  * \code{rank} the matrix rank of \code{object} or its returned inverse
#'  (\code{rank = 1} for any vector).
#'
#' @seealso \link[bigpcor]{inverse.matrix} for the default method.
#'
inverse <- function (object, ...) {
  inverse.matrix (X = object, ...)
}


# Define a generic function inverse
setGeneric(name = "inverse",
           def = inverse)
