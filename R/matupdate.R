# Deprecated code
# DELETE THE FILE

#
# Updating Objects Resulting From Matrix Operations
#
# \code{matupdate} is a generic function whose default method updates a square
# matrix, inverse of an inner product matrix \eqn{t(X) \times X}, after
# the \eqn{n \times p} matrix \eqn{X} has been updated either by including an
# additional column \eqn{v} into \eqn{X}, or by removing one column from
# \eqn{X}.
#
# @param object a numeric matrix-like object of a class with a
# \code{matupdate} method. For the default method, \code{object} is a numeric
# matrix \eqn{A} of size \eqn{p \times p}, inverse of the inner product matrix
# \eqn{t(X) \times X}. \code{object} is of the class \code{"inner.inv"} which
# is a matrix with a named attribute 'X' (that can be \code{NULL}, i.e. missing).
#
# @param add logical, adding information into \eqn{X}?
#
# @param position a positive integer, indicating where exactly is the
# information being added or dropped. For the default method, indicates the
# position \eqn{j} of:
# \itemize{
#  \item{the vector \eqn{v} being added into \eqn{X} when \code{add = TRUE}}{, here \eqn{1 \leq j \leq p + 1} with default value \code{position}=\eqn{j}=\code{p+1};}
#  \item{the vector \eqn{v} being dropped from \eqn{X} when \code{add = FALSE}}{, here \eqn{1 \leq j \leq p} with default value \code{position}=\eqn{j}=\code{p}.}
# }
#
# @param v a numeric matrix-like object providing new information to be used
# to update \code{object}. For the default method, a length \eqn{n} numeric
# vector, required if \code{add = TRUE} (the default) and not required nor
# used when \code{add = FALSE}.
#
# @param X a numeric matrix-like object, contains the basic information used
# to obtain \code{object}. For the default method, a \eqn{n \times p} numeric
# matrix.
#
# @param ... additional argument passed to or from other methods.
# Currently not used by the default method.
#
# @export matupdate
# @exportMethod matupdate
#
# @details
# Additional \code{matupdate} methods can be defined for specific classes.
#
# @return For the default method \link{matupdate.inner.inv}, an object of class
# \code{"inner.inv"}, i.e. a matrix with a named attribute \code{"X"} (which
# can be \code{NULL}, i.e. missing).
#
# @seealso \link[bigpcor]{matupdate.inner.inv}.
#
#matupdate <- function(object, add = TRUE, position = NULL, v, X, ...) {
#  matupdate.inner.inv (object = object, add = add,
#                     position = position,
#                      v = v, X = X, ...)
#}


# Define a generic function matupdate
#setGeneric(name = "matupdate",
#           def = matupdate)
