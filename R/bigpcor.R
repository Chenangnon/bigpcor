# The package uses a relatively fast approach for calculations,
# taking advantage linear algebra and library 'bigcor'.
#
#' Pairwise Partial Correlations
#'
#' Calculate partial correlations between elements of two sets of variables
#' (\code{x} and \code{y}) conditioning on another set of variables (\code{z}),
#' e.g. confounders to control for.
#'
#' @param x,y integer vector, or vector of labels (column names of \code{data} as
#' returned by \link[base]{colnames}(\code{data})), or matrix or data frame giving
#' variables whose partial correlations are desired.
#'
#' Unlike for \link[stats]{cor}, the default for \code{y} (i.e. \code{NULL}) is not
#' generally equivalent to \code{y = x}. Here, \code{y = NULL} is equivalent to
#' \code{y = x} only when \code{data} is not \code{NULL} (see section "Details").
#'
#' When \code{data} is specified and not \code{NULL}, \code{x} must be an integer
#' vector or a vector of labels; and \code{y} must also be one of an integer vector, a
#' vector of labels, or \code{NULL}; both \code{x} and \code{y}, when not
#' \code{NULL}, giving column indices/labels indicating variables in \code{data}.
#'
#' @param z integer vector, vector of labels (column names of \code{data}),
#' or matrix or data frame. \code{z} specifies the conditioning set for
#' computing partial correlations. When \code{data} is specified and not
#' \code{NULL}, \code{z} must be an integer vector, or a vector of labels,
#' or \code{NULL}.
#'
#' Note that \code{z = NULL} does not correspond to an empty conditioning set,
#' but rather conditioning on all other variables in \code{x} and \code{y}, i.e.
#' \code{z = NULL} is equivalent to \code{z = c(x, y)} when \code{data} is not
#' \code{NULL}. An empty conditioning set corresponds to marginal correlation,
#' consider the \R routine \link[stats]{cor} for that purpose, or use, for
#' instance, the function \link[propagate]{bigcor} of library \code{propagate}
#' for very large datasets.
#'
#' @param data matrix or data frame. When specified and not \code{NULL},
#' \code{x}, \code{y} and \code{z} are interpreted as column indices or labels
#' (column names of \code{data}) of corresponding numeric variables in \code{data}.
#' In this case, \code{x} must be an integer vector or a vector of labels, and
#' each of \code{y} and \code{z} must also be an integer vector, vector of
#' labels, or \code{NULL}.
#'
#' @param use an optional character string giving a method for computing
#' correlations in the presence of missing values. This must be (an abbreviation
#' of) one of the strings \code{"everything"}, \code{"all.obs"},
#' \code{"complete.obs"}, \code{"na.or.complete"}, or \code{"pairwise.complete.obs"}.
#'
#' @param method a character string indicating which correlation coefficient
#' is to be computed. One of \code{"pearson"} (default), \code{"kendall"}, or
#' \code{"spearman"}: can be abbreviated.
#'
#' @param blocksize the size of sub-matrices for block-wise calculation of
#' large covariance matrices (see argument \code{size} of the function
#' \link[propagate]{bigcor} of library \code{propagate}). Defaults to 2000.
#'
#' @param cl a cluster object, created by package \code{parallel} or package
#' \code{snow}. Defaults to the registered default cluster obtained as
#' \code{cl = parallel::getDefaultCluster()}. Only used when \code{z} is not
#' \code{NULL}.
#'
#' @param chunk.size scalar number, number of units for scheduling parallel
#' tasks. Only used when \code{cl} is a proper cluster object.
#'
#' @param verbose logical. If \code{TRUE}, information is printed in the
#' console when running on a large data matrix (100 variables or more).
#'
#' @param ... additional argument passed to or from other methods. Currently not
#' used.
#'
#' @details
#' Partial correlations express the specific portion of variance explained after
#' controlling for the effect of possible confounding variables when assessing
#' the correlation between two variables \insertCite{kim2015ppcor}{bigpcor}.
#'
#' The function takes advantage of:
#' \itemize{
#'  \item{}{function \link[propagate]{bigcor} from library \code{propagate}
#'  \insertCite{Spiess2018propagate}{bigpcor} to compute covariance matrix
#'  for large datasets (100 variables or more);}
#'  \item{}{linear algebra (\insertCite{khan2008updating;textual}{bigpcor}'s
#'  formula for rank-1 update of inverse of inner product) to efficiently
#'  compute partial correlations by reducing repeated matrix inverse
#'  computations to one inversion and then simpler matrix and vector
#'  multiplications.}
#' }
#'
#' The conditioning set considered to compute partial correlations depends on
#' the specified arguments \code{x}, \code{y}, and \code{z}.
#'
#' \itemize{
#'  \item{\code{x} and \code{z} are specified: }{\code{bigpcor} returns partial
#'  correlations \code{cor(x_i,y_j)} given all variables (columns) in \code{z}
#'  (or the corresponding columns in \code{data}). If \code{y = NULL},
#'  \code{y = x} is used. Specifying \code{z = c(x, y)} is equivalent to
#'  \code{z = NULL} only when \code{data} is not \code{NULL} and \code{x} and
#'  \code{y} are integer vectors or vectors of labels.}
#'  \item{Only \code{x} and \code{y} are specified: }{\code{bigpcor} returns
#'  partial correlations \code{cor(x_i,y_j)} given all other variables (columns)
#'  in \code{x} and \code{y}. Specifying \code{y = x} is equivalent to
#'  \code{y = NULL} only when \code{data} is not \code{NULL} and \code{x} and
#'  \code{y} are integer vectors or vectors of labels (indeed, duplicates
#'  appearing in both the pair and the conditioning set can be removed in this
#'  case). In the general case, since only pairs with one element in \code{x}
#'  and the other in \code{y} are considered, appropriately specifying \code{x}
#'  and \code{y} can help avoid computing all pairwise partial correlations
#'  between all variables (in \code{x} and \code{y}) if only a subset is desired.}
#'  \item{Only \code{x} is specified: }{\code{bigpcor} returns partial correlations
#'   \code{cor(x_i,x_j)} between columns of \code{x} (\code{x_i} and \code{x_j}
#'   are respectively the \code{i}th and \code{j}th columns of \code{x}, or the
#'   corresponding columns in \code{data}) given all other variables (columns)
#'   in \code{x}.}
#' }
#'
#' When \code{x,y,z} are column indices (numeric vectors), they must be
#' \code{"integers"} (checked using \code{floor(x) == x}). When \code{x,y,z}
#' are column indices or column labels of \code{data} (column names \code{data}
#' as returned by \link[base]{colnames}), duplicates are removed from each set
#' using the function \link{unique}. Thus, to avoid unnecessarily dealing with
#' singular matrix inversion, we recommend to specify the argument \code{data}
#' and give \code{x}, \code{y} and \code{z} as column indices or labels.
#'
#' @return A matrix of partial correlations, rows corresponding to \code{x} and
#' columns to \code{y} (replaced by \code{x} if \code{NULL}). The
#' returned matrix has two named attributed:
#'   * \code{S}: the size of the used conditioning set;
#'   * \code{rank}: the matrix rank of the used conditioning set.
#'
#' @export bigpcor
#' @importFrom propagate bigcor
#' @importFrom stats cor
#' @importFrom stats cov
#' @importFrom stats cov2cor
#' @importFrom utils combn
#' @importFrom stats na.omit
#' @importFrom parallel getDefaultCluster
#' @importFrom Rdpack insert_all_ref
#' @importFrom Rdpack insert_ref
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \link[ppcor]{spcor} of library \code{ppcor} for semi partial correlations,
#' \link[stats]{cor} for marginal correlation,
#' and \link[propagate]{bigcor} of library \code{propagate} for marginal
#' correlations from large datasets.
#'
#' @examples
#' # This example mimicks the selection of confounders of gene expressions
#' # in a genomic network. We consider a simulated data with 50 genes ('fiftygenes').
#'
#' #* Task: compute pairwise partial correlations between expressions on one hand
#' # and candidate confouders on the other hand, conditioning on genetic variants.
#'
#' #* Set number of genetic variants (q), genes (p) and confounders (u) in the
#' # dataset 'fiftygenes'
#' p <- q <- 50
#' u <- 100
#'
#' #* Using 'bigpcor'
#' Time <- system.time({
#'   res <- bigpcor (x = fiftygenes[, (q+1):(q+p)],
#'                   y = fiftygenes[, (q+p+1):(q+p+u)],
#'                   z = fiftygenes[, 1:q])
#' })
#'
#' # A glance at the results for the first five genes and five confounders
#' res[1:5, 1:5]
#'
#' # Elapsed time (seconds)
#' Time
#'
#' # Library 'ppcor' is required to run the remaining of the example
#'
#' \donttest{
#' ##* Comparing 'bigpcor' with a repeated call to 'pcor' of library "ppcor"
#' # Install 'ppcor' if not installed
#' if (!"ppcor" %in% rownames(installed.packages()))
#'   install.packages('ppcor', dependencies = TRUE)
#' require(ppcor)
#'
#' # Define a function to disable the display of warnings from 'ppcor::pcor'
#' # This is Martin Maechler's 'tryCatch.W.E' which catches and saves both errors
#' # and warnings, and in the case of a warning, also keeps the computed result.
#' catch.conditions <- function (expr) {
#'   W <- NULL
#'   w.handler <- function(w) {
#'     W <<- w
#'     invokeRestart("muffleWarning")
#'   }
#'   list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
#'                                   warning = w.handler), warning = W)
#' }
#'
#' #* Use the 'pcor' function of library "ppcor"
#' # The following took about 5 seconds on a MAC OS Ventura 13.4.1 system with 16 GB RAM
#' W.E. <- catch.conditions({
#'   Time0 <- system.time({
#'     res0 <- apply(expand.grid(E = (q+1):(q+p), C = (q+p+1):(q+p+u)),# Combinations of E and C
#'                   MARGIN = 1,
#'                   FUN = function(x) { # call 'ppcor::pcor'
#'                     pcor (fiftygenes[, c(x, 1:q)])$estimate[1,2]
#'                   })
#'   })
#' })
#'
#' # Run the following to display catched warning(s)
#' # W.E.$warning # Not run
#'
#' # Organize res0 into a matrix (rows for genes and columns for confounders)
#' res0 <- matrix(res0, nrow = p, ncol = u, byrow = FALSE)
#'
#' # A glance at the results for the first five genes and five confounders
#' res0[1:5, 1:5]
#'
#' # Compare results to 10 decimal places
#' sum(round(res, 10) != round(res0, 10)) # zero means no difference
#'
# Compare elapsed times (seconds)
#' Time0
#' Time0[3]/Time[3] # 'bigpcor' is about four fold faster than 'ppcor' in this example
#' }
#'
#\donttest{
# # Set number of variants (q), genes (p) and confounders in the dataset 'hundredgenes'
# p <- q <- 100
# u <- 200
#
# #* Using 'bigpcor'
# Time <- system.time({
#   res <- bigpcor (x = hundredgenes[, (q+1):(q+p)],
#                   y = hundredgenes[, (q+p+1):(q+p+u)],
#                   z = hundredgenes[, 1:q])
# })
#
# # A glance at the results for the first five genes and five confounders
# res[1:5, 1:5]
#
# # The following took about 150 seconds on a MAC OS Ventura 13.4.1 system with 16 GB RAM
# W.E. <- catch.conditions({
#   Time0 <- system.time({
#     res0 <- apply(expand.grid(E = (q+1):(q+p), C = (q+p+1):(q+p+u)),# Combinations of E and C
#                   MARGIN = 1,
#                   FUN = function(x) { # call 'ppcor::pcor'
#                     pcor (hundredgenes[, c(x, 1:q)])$estimate[1,2]
#                   })
#   })
# })
#
# # Run the following to display catched warning(s)
# # W.E.$warning # Not run
#
# # Organize res0 into a matrix (rows for genes and columns for confounders)
# res0 <- matrix(res0, nrow = p, ncol = u, byrow = FALSE)
#
# # A glance at the results for the first five genes and five confounders
# res0[1:5, 1:5]
#
# # Compare results to 10 decimal places
# sum(round(res, 10) != round(res0, 10)) # zero means no difference
#
# # Compare elapsed times
# Time0
# Time0[3]/Time[3] # 'bigpcor' is about ten fold faster than 'ppcor' in this example
# }
#'
bigpcor <- function (x, y = NULL, z = NULL, data = NULL,
                     use = "everything",
                     method = c("pearson", "spearman", "kendall"),
                     blocksize = NULL,
                     cl = parallel::getDefaultCluster(),
                     chunk.size = NULL, verbose = FALSE, ...) {
  ### Handle arguments 'method' and 'use'
  ## use
  na.method <- pmatch(use[1], c("all.obs", "complete.obs", "pairwise.complete.obs",
                                "everything", "na.or.complete"))
  if (is.na(na.method))
    stop("invalid 'use' argument")

  ## method
  method <- match.arg(method)

  ### Ensure symmetry of the function with respect to x and y
  if (missing(x) | is.null(x)) {
    if (is.null(y))
      stop("at least one of 'x' and 'y' must be specified and non NULL")
    x <- y
    y <- NULL
  }

  ### Check and organize x, y, z, and data argument
  p <- q <- r <- dnames <- NULL
  eval(check.pcor.data())

  ### Case the conditioning set is NULL, condition on all other variables
  if (r == 0) {
    ### Block size if required
    if (p + q > 99) {
      if (is.null(blocksize)) {
        blocksize <- min(p+q, 2000)
      }
      else {
        blocksize <- min(p+q, max(blocksize, 10))
      }
    }

    if (is.null(y)) {
      # Get the covariance matrix
      if (p > 99) {
        if (verbose) {
          covxy <- propagate::bigcor(x = data[,x], fun = "cov",
                                     size = blocksize, verbose = verbose,
                                     use = use, method = method)
        }
        else {
          covxy <- catch.conditions({
            propagate::bigcor(x = data[,x], fun = "cov",
                              size = blocksize, verbose = verbose,
                              use = use, method = method)
          })$value
        }
        covxy <- covxy[1:p, 1:p]
      }
      else {
        covxy <- cov(x = data[,x], use = use, method = method)
      }

      # Transform into partial correlations
      H <- mp.inverse(covxy)
      Sigmarank <- attr(H, "rank")
      H <- - cov2cor(H)
      diag(H) <- 1
      if (is.null(dnames))
        dnames <- paste0('x', x)
      dimnames(H) <- list(dnames, dnames)

      # Return with the size of the conditioning set as an attribute
      out <- structure(H, S = p - 2, rank = Sigmarank - 2)

      return(out)
    }
    else {
      # Get the covariance matrix
      xy <- unique(c(x,y))
      nxy <- length(xy)
      if (p + q > 99) {
        if (verbose) {
          covxy <- propagate::bigcor(x = data[, xy, drop = FALSE], fun = "cov",
                                     size = blocksize, verbose = verbose,
                                     use = use, method = method)
        }
        else {
          covxy <- catch.conditions({
            propagate::bigcor(x = data[, xy, drop = FALSE], fun = "cov",
                              size = blocksize, verbose = verbose,
                              use = use, method = method)
          })$value
        }
        covxy <- covxy[1:nxy, 1:nxy]
      }
      else {
        covxy <- cov(x = data[, xy, drop = FALSE], use = use, method = method)
      }

      # Transform into partial correlations
      H <- mp.inverse(covxy)
      Sigmarank <- attr(H, "rank")
      H <- - cov2cor(H)
      if (is.null(dnames))
        dnames <- c(paste0('x', x),
                    paste0('y', y))
      dimnames(H) <- list(dnames, dnames)

      # Return only desired rows and columns
      out <- structure(H[x, y], S = nxy - 2, rank = Sigmarank - 2)

      return(out)
    }
  }

  ### Cases where we have variables in 'z' to condition on

  ## Center data columns to be used
  data[, c(x,y,z)] <- scale(data[, c(x,y,z), drop = FALSE],
                            center = TRUE, scale = FALSE)

  ## Inverse covariance matrix of the conditioning set
  ## (any required data transformation is taken care of by 'stats::cov' or 'propagate::bigcor')
  if (p > 99) {  # If 100 variables or more in the conditioning set
    # Block size
    if (is.null(blocksize)) {
      blocksize <- min(r, 2000)
    }
    else {
      blocksize <- min(r, max(blocksize, 10))
    }

    # Compute covariance matrix
    if (verbose) {
      invcovz <- propagate::bigcor(x = data[, z, drop = FALSE], fun = "cov",
                                   size = blocksize, verbose = verbose,
                                   use = use, method = method)
    }
    else {
      invcovz <- catch.conditions({
        propagate::bigcor(x = data[, z, drop = FALSE], fun = "cov",
                          size = blocksize, verbose = verbose,
                          use = use, method = method)
      })$value
    }
    invcovz <- invcovz[1:r, 1:r]

    # Compute inverse covariance
    invcovz <- mp.inverse(invcovz)
    Sigmarank <- attr(invcovz, "rank")
    FullRank <- attr(invcovz, "m") - Sigmarank == 0
  }
  else { # For 99 or less variables in the conditioning set
    invcovz <- mp.inverse(cov(x = data[, z, drop = FALSE], use = use, method = method))
    Sigmarank <- attr(invcovz, "rank")
    FullRank <- attr(invcovz, "m") - Sigmarank == 0
  }

  ## rescale 'invcovz' (we only need 't(X) %*% X', not the division by n - 1 as done by 'cov')
  invcovz <- invcovz / (NROW(data) - 1) # Divide because we inverted 'cov(data[,z])', so 'invcovz' is what we want times (n - 1)

  ## Data transform if required  (these steps are after the function 'stats::cor')
  if (method[1] %in% c("spearman", "kendall")) {
    if (na.method %in% c(2L, 5L)) {
      nas <- attr(na.omit(data), "na.action")
      data <- Rrank (RdropNA(data, nas))
    }
    else if (na.method != 3L) {
      data <- Rrank (data)
    }
    else stop("cannot handle 'pairwise.complete.obs'")
  }
  else {
    if (na.method %in% c(2L, 5L)) {
      nas <- attr(na.omit(data), "na.action")
      data <- RdropNA(data, nas)
    }
    else if (na.method == 3L)
      stop("cannot handle 'pairwise.complete.obs'")
  }

  ## Handle case of missing 'y' (use y = x;
  #                              but only compute lower triangle elements
  #                              of the results, and use symmetry to fill
  #                              the upper triangle)
  if (is.null(y)) {
    ## p must be 2 or more
    stopifnot(p > 1)

    ## Create a grid of (x,y) pairs (i.e. pairs of indices of columns in data)
    ## Here, we only consider indices (i, j) such that j < i
    xygrid <- t(combn(x, m = 2))

    ## Get the partial correlations (lower triangle elements)
    pcorrs <- matteApply (X = cbind(xygrid, FALSE),
                          MARGIN = 1,
                          FUN = pcorcxyz,
                          invcovz = invcovz,
                          z = z,
                          fullzrank = FullRank,
                          data = data,
                          chunk.size = chunk.size, cl = cl)

    ## Reorganize into a matrix of partial correlations
    H <- diag(p)
    H[lower.tri(H, diag = FALSE)] <- pcorrs # (lower triangle)

    ## Make H symmetric
    H <- H + t(H) # (upper triangle)
    diag(H) <- 1

    ## Name result rows and columns
    if (!is.null(dnames))
      dimnames(H) <- list(dnames[x], dnames[x])
    else {
      dimnames(H) <- list(paste0('x', x), paste0('x', x))
    }

    out <- structure(H, S = r, rank = Sigmarank)

    return(out)
  }

  ## Create a grid of (x,y) pairs (i.e. pairs of indices of columns in data)
  xygrid <- as.matrix(expand.grid(x = x, y = y))
  grid <- 1:nrow(xygrid) # Position in the grid
  pcor_position <- grid  # Position of the partial correlation

  ## Remove duplicates if any
  intxy <- intersect(x, y)
  if (length(intxy)) {
    intxy <- (xygrid[,1] %in% intxy) + (xygrid[,2] %in% intxy) == 2
    if (any(intxy)) {
      pcor_position[intxy] <- matteSapply (X = pcor_position[intxy],
                                           FUN = remove.duplicates,
                                           p = p, xygrid = xygrid,
                                           chunk.size = chunk.size, cl = cl)
    }
  }

  ## Mark duplicated xygrid values
  gduplicates <- grid != pcor_position

  ## Get the partial correlations
  pcorrs <- matteApply (X = cbind(xygrid, gduplicates),
                        MARGIN = 1,
                        FUN = pcorcxyz,
                        invcovz = invcovz,
                        z = z,
                        fullzrank = FullRank,
                        data = data,
                        chunk.size = chunk.size, cl = cl)

  ## Reorganize into a matrix of partial correlations
  H <- matrix(0, nrow = p, ncol = q)
  xygrid0 <- as.matrix(expand.grid(x = 1:p, y = 1:q))
  H[xygrid0] <- pcorrs[pcor_position]

  ## Name result rows and columns
  if (!is.null(dnames))
    dimnames(H) <- list(dnames[x], dnames[y])
  else {
    dimnames(H) <- list(paste0('x', x), paste0('y', y))
  }

  out <- structure(H, S = r, rank = Sigmarank)

  return(out)
}

check.pcor.data <- function(x, y, z, data) {
  expression({
    if (is.null(data)) {
      # x
      x <- as.matrix(x)
      stopifnot(is.numeric(x))
      p <- ncol(x)
      n <- NROW(x)

      # y
      if (is.null(y)) {
        q <- 0
      }
      else {
        y <- as.matrix(y)
        stopifnot(is.numeric(y))
        stopifnot(NROW(y) == n)
        q <- ncol(y)
      }

      # z
      if (is.null(z)) {
        # Number of variables in z
        r <- 0
      }
      else {
        z <- as.matrix(z)
        stopifnot(is.numeric(z))
        stopifnot(NROW(z) == n)

        # Size of the conditioning set
        r <- ncol(z)
      }

      # bind the three sets
      data <- cbind(x, y, z)
      m <- p + q + r
      xnames <- colnames(x)
      ynames <- if (q > 0) colnames(y)
      znames <- if (r > 0) colnames(z)
      x <- 1:p
      y <- if (!is.null(y)) (p+1):(p+q)
      z <- if(!is.null(z)) (p+q+1):(p+q+r)

      if (is.null(xnames))
        xnames <- paste0('x', x)
      if (!is.null(y) & is.null(ynames))
        ynames <- paste0('y', y)
      if (!is.null(z) & is.null(znames))
        znames <- paste0('z', z)
      dnames <- c(xnames, ynames, znames)
    }
    else {
      # data
      stopifnot(!is.null(dim(data)))
      m <- ncol(data)
      n <- NROW(x)
      dnames <- colnames(data)

      # x
      stopifnot(is.vector(x), (is.numeric(x) | is.character(x)))
      if (is.numeric(x)) {
        stopifnot(all(x == floor(x)))
        stopifnot(all(c(x > 0, x <= m)))
        x <- unique(x)
      }
      else {
        if (is.null(dnames)) {
          stop("'data' must have non null colnames when 'x' contains labels (strings)")
        }
        stopifnot(all(x %in% dnames))
        x <- unique(x)
        x <- sapply(x, FUN = function(i) which(dnames == i))
      }
      p <- length(x)

      # y
      if (is.null(y)) {
        q <- 0
      }
      else {
        stopifnot(is.vector(y), (is.numeric(y) | is.character(y)))
        if (is.numeric(y)) {
          stopifnot(all(y == floor(y)))
          stopifnot(any(c(y > 0, y <= m)))
          y <- unique(y)
        }
        else {
          if (is.null(dnames)) {
            stop("'data' must have non null colnames when 'y' contains labels (strings)")
          }
          stopifnot(all(y %in% dnames))
          y <- unique(y)
          y <- sapply(y, FUN = function(i) which(dnames == i))
        }
        q <- length(y)
      }

      # z
      if (is.null(z)) {
        # Number of variables in z
        r <- 0
      }
      else {
        stopifnot(is.vector(z), (is.numeric(z) | is.character(z)))
        if (is.numeric(z)) {
          stopifnot(all(z == floor(z)))
          stopifnot(any(c(z > 0, z <= m)))
          z <- unique(z)
        }
        else {
          if (is.null(dnames)) {
            stop("'data' must have non null colnames when 'z' contains labels (strings)")
          }
          stopifnot(all(z %in% dnames))
          z <- unique(z)
          z <- sapply(z, FUN = function(i) which(dnames == i))
        }

        # Size of the conditioning set
        r <- length(z)
      }

      # Ensure all indexed columns in data are numerical
      stopifnot(is.numeric(as.matrix(data[, c(x,y,z), drop = FALSE])))

      # Get data column names
      if (is.null(dnames)) {
        dnames <- paste0('V', 1:m)
        colnames(data) <- dnames
      }
      else {
        dnames <- colnames(data) <- make.unique(dnames)
      }
    }
  })
}

## Child function: workhorse for 'bigpcor'
pcorcxyz <- function (ij, invcovz, z, data, fullzrank = TRUE) {
  ## return NA a duplicate
  if (is.na(ij[3])) {
    return(NA)
  }

  ## correlation is one if xi and yj are rhe same
  ij <- ij[1:2]
  if (ij[1] == ij[2]) {
    return(1)
  }

  dup <- which(z %in% ij)
  ## If xi and yj are already in z
  if (length(dup) == 2) { # elements in z are unique in z, so length(dup) == 2 implies that both xi and yj are in z
    H <- - cov2cor(invcovz)
    #diag(H) <- 1
    return(H[dup[1], dup[2]])
  }

  ## If one of xi and yj is already in z
  if (length(dup) == 1) {
    invcovz <- add.last.inner.inv(invcovz, X = data[,z, drop = FALSE], v = data[, ij[!(ij %in% z)]])
    H <- - cov2cor(invcovz)
    #diag(H) <- 1
    return(H[dup, length(z) + 1])
  }

  ## If none of xi and yj is already in z
  # First include i into invcovz
  invcovz <- add.last.inner.inv(invcovz, X = data[,z, drop = FALSE],
                                v = data[, ij[1]], MP = !fullzrank)

  # Then include j
  invcovz <- add.last.inner.inv(invcovz, X = data[,c(z, ij[1]), drop = FALSE],
                                v = data[, ij[2]], MP = !fullzrank)
  H <- - cov2cor(invcovz)
  #diag(H) <- 1
  return(H[length(z) + 1, length(z) + 2])

}

## Child to remove duplicates in 'x' and 'y' when 'z' is not NULL
remove.duplicates <- function (j, p, xygrid) {
  if (j <= p)
    return(j) # Because no duplicate can happen for the combination of all elements of x with the first element of y
  xeq <- xygrid[j,1] == xygrid[1:(j-1), 2]
  yeq <- xygrid[j,2] == xygrid[1:(j-1), 1]
  xyeq <- which((xeq + yeq) == 2)
  if (length(xyeq)) {
    return(xyeq[1])
  }
  else {
    return(j)
  }
}

## Child from function 'cor' of library 'stats'
Rrank <- function(u) {
  if (length(u) == 0L)
    u
  else if (is.matrix(u)) {
    if (nrow(u) > 1L)
      apply(u, 2L, rank, na.last = "keep")
    else row(u)
  }
  else rank(u, na.last = "keep")
}

## Child from function 'cor' of library 'stats'
RdropNA <- function(x, nas) {
  if (length(nas)) {
    if (is.matrix(x))
      x[-nas, , drop = FALSE]
    else x[-nas]
  }
  else x
}
