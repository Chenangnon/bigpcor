#' Simulated Fake Genomic Network Data
#'
#' Datasets mimicking genomic networks.
#' This is aimed at showcasting the selection of confounders of gene expressions
#' using package \code{bigpcor}.
#'
#' @details
#' Four simulated datasets are available.
#' Each dataset has three sets of variables,
#'
#'  * a set 'V' with \code{p} variables, representing bi-allelic genetic variants
#'  * a set 'E' with \code{p} variables, representing gene expression values
#'  * a set 'C' with \code{u} variables, representing candidate confounders
#'                                       for the expression values in 'E'.
#'
#'  For \code{tengenes}, \code{p=10} and \code{u = 20} with sample size \eqn{n=50}.
#'
#'  For \code{fiftygenes}, \code{p=50} and \code{u = 100} with sample size \eqn{n=300}.
#'
#'  For \code{hundredgenes}, \code{p=100} and \code{u = 200} with sample size \eqn{n=500}.
#'
# ## Sampling the three set of variables
# #Numbers of variables and sample size for each sets
# q <- 100
# p <- 100
# u <- 200
# n <- 500
#
# # set seed for reproducibility
# set.seed(100)
#
# # For simplicity, all variables are mutually independent (the sampled expressions
# #     are independent of each other and independent of the q genetic variants)
# #     'V' taken from a discrete uniform over {0,1,2},
# #     'E' and 'C' from standard normal distribution.
#
# # The q variants
# V <- matrix(sample(0:2, size = q * n, replace = TRUE),
#                    nrow = n,
#                    ncol = q)
# colnames(V) <- paste0("V", 1:q)
#
# # The p expressions
# E <- matrix(rnorm(p * n),
#             nrow = n,
#             ncol = p)
# colnames(E) <- paste0("E", 1:p)
#
# # The u confounders
# C <- matrix(rnorm(u * n),
#                   nrow = n,
#                   ncol = u)
# colnames(C) <- paste0("C", 1:u)
#
#
#' @docType data
#'
#' @aliases tengenes
#' @aliases fiftygenes
#' @aliases hundredgenes
#'
#' @usage data(tengenes)
#'
#' @usage data(fiftygenes)
#'
#' @usage data(hundredgenes)
#'
#' @format An object of class \code{"data.frame"}.
#'
#' @keywords datasets
#'
#' @source Simulated
#'
#' @examples
#' data(tengenes)
#' \donttest{
#' pairs(tengenes[,11:20])
#' }
"tengenes"
"fiftygenes"
"hundredgenes"
