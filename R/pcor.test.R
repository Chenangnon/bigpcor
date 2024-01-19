#
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats p.adjust
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom stats pt
# @importFrom stats cov
#' @importFrom stats cor
#
# KEEP 'qvalue' option only after writing a self-qvalue function
## As package pvalue is no longer on CRAN, and I do not want DIRECT dependence on Bioconductor maintained packages
## OR find an ALTERNATIVE TO qvalue correction (Addis, ...)
pcor.test <- function (x, y = NULL, z = NULL,
                       data = NULL, subset, na.action = na.omit,
                       method = c("pearson", "kendall", "spearman"),
                       blocksize = NULL,
                       rho = 0, statistic = c("t", "z"),
                       alternative = c("two.sided", "less", "greater"),
                       p.adjust.method = "none", pi0.method = "smoother",
                       fdr.level = 1 - conf.level,
                       exact = NULL, conf.level = 0.95, continuity = FALSE,
                       cl = NULL, chunk.size = NULL, verbose = FALSE, ...) {
  ### Control arguments
  # method
  method <- match.arg(method)

  # alternative
  alternative <- match.arg(alternative)

  # p.adjust.method
  p.adjust.method <- match.arg(p.adjust.method,
                               choices = c("none", "holm", "hochberg", "hommel",
                                           "bonferroni", "BH", "BY", "fdr",
                                           "qvalue"))

  # rho
  stopifnot(is.null(rho), all(abs(rho) < 1))

  # Test statistic
  if (missing(statistic)) {
    if (all(rho == 0)) {
      # use t-test as in 'cor.test' if the null hypothesis is: no correlation
      statistic <- "t"
    }
    else {
      # Otherwise, use z-test (Fisher's transform) as in 'cor.test'
      statistic <- "z"
    }
  }
  else {
    statistic <- match.arg(statistic)
  }

  ### Names of the inputs in the calling environment
  DNAME <- paste(if (!is.null(x)) deparse1(substitute(x)), if (is.null(z)) "and" else ",", if (!is.null(y)) deparse1(substitute(y)), if (!is.null(z)) "and", deparse1(substitute(z)))

  ## Check and organize x, y, z, and data argument
  eval(check.pcor.data())

  ### Handle arguments 'subset' and 'na.action' if supplied
  data <- model.frame(formula = ~ .,
                      data = as.data.frame(data),
                      subset = subset, na.action = na.action)[, -1, drop = FALSE]
  n <- NROW(data)

  ### Use bigpcor to obtain the partial correlation estimates
  ESTIMATE <- bigpcor(x = x, y = y, z = z, data = data,
                      method = method, blocksize = blocksize,
                      cl = cl, chunk.size = chunk.size,
                      verbose = verbose, ...)

  ### Size of the conditioning set
  S <- attr(ESTIMATE, "S")

  ### Degrees of Freedom for test
  df0 <- n - S - 3L
  switch(statistic,
         t = {
           TTEST <- TRUE
           df <- df0
         },
         {
           TTEST <- FALSE
           df <- Inf
         })

  ### Compute statistics
  NVAL <- 0
  conf.int <- FALSE
  if (method == "pearson") {
    if (df0 < 0)
      stop("not enough finite observations")
    method <- "Pearson's product-moment correlation"
    names(NVAL) <- "correlation"
    PARAMETER <- df
    # Standard error of estimates
    sigma <- 1/sqrt(df0)

    if (TTEST) {
      STATISTIC <- sqrt(df0) * ESTIMATE/sqrt(1 - ESTIMATE^2)
      attr(STATISTIC, "statistic") <- "t"

      # Calculate p-values
      PVAL <- switch(alternative,
                     less = pt(STATISTIC, df = df, lower.tail = TRUE),
                     greater = pt(STATISTIC, df = df, lower.tail = FALSE),
                     two.sided = 2 * pt(abs(STATISTIC), df = df, lower.tail = FALSE))

      # Adjust p-values if required
      switch(p.adjust.method,
             qvalue = { # If q-value correction
               PVAL <- matrix(q.adjust(c(PVAL),
                                       fdr.level = NULL,
                                       pfdr = FALSE,
                                       lfdr.out = TRUE,
                                       pi0.method = pi0.method),
                              nrow = p, ncol = q, byrow = FALSE)
             },
             { # If Multiple Comparisons Methods
               PVAL <- matrix(p.adjust(c(PVAL),
                                       method = p.adjust.method,
                                       n = length(PVAL)),
                              nrow = p, ncol = q, byrow = FALSE)
             }
             )

      # Build confidence interval
      if (df0 > 0) {
        if (!missing(conf.level)) {
          if(length(conf.level) != 1 || !is.finite(conf.level) || conf.level < 0 || conf.level > 1)
            stop("'conf.level' must be a single number between 0 and 1")
        }
        conf.int <- TRUE

        # Fishers' z-transform
        fz <- atanh(ESTIMATE)

        cint <- switch(alternative,
                       two.sided = sigma * qnorm((1 + conf.level)/2),
                       sigma * qnorm(conf.level))

        cint <- switch(alternative,
                       less = list(lower = -1, upper = tanh(fz + cint)),
                       greater = list(lower = tanh(fz - cint), upper = 1),
                       two.sided = list(lower = tanh(fz - cint),
                                        upper = tanh(fz + cint)))
        attr(cint, "conf.level") <- conf.level
      }
    }
    else {
      fz <- atanh(ESTIMATE)
      STATISTIC <- fz - atanh(rho)
      attr(STATISTIC, "statistic") <- "z"

      if (df0 > 0) {
        if (!missing(conf.level)) {
          if(length(conf.level) != 1 || !is.finite(conf.level) || conf.level < 0 || conf.level > 1)
            stop("'conf.level' must be a single number between 0 and 1")
        }
        conf.int <- TRUE



        cint <- switch(alternative,
                       two.sided = sigma * qnorm((1 + conf.level)/2),
                       sigma * qnorm(conf.level))

        cint <- switch(alternative,
                       less = list(lower = -1, upper = tanh(fz + cint)),
                       greater = list(lower = tanh(fz - cint), upper = 1),
                       two.sided = list(lower = tanh(fz - cint),
                                        upper = tanh(fz + cint)))
        attr(cint, "conf.level") <- conf.level
      }

      # Calculate p-values
      PVAL <- switch(alternative,
                     less = pnorm(STATISTIC, lower.tail = TRUE),
                     greater = pnorm(STATISTIC, lower.tail = FALSE),
                     two.sided = 2 * pnorm(abs(STATISTIC), lower.tail = FALSE))

    }
  }
  else {
    if (n < 2)
      stop("not enough finite observations")
    PARAMETER <- NULL

    TIES <- outer(apply(data[, x, drop = FALSE], MARGIN = 2, FUN = function(x) length(unique(x))),
                  apply(data[, y, drop = FALSE], MARGIN = 2, FUN = function(y) length(unique(y))),
                  FUN = function(x,y) min(x, y)) < n

    if (method == "kendall") {
      method <- "Kendall's rank correlation tau"
      names(NVAL) <- "tau"
      NA.EST <- !is.finite(ESTIMATE)
      if (any(NA.EST)) {
        ESTIMATE[NA.EST] <- STATISTIC[NA.EST] <- PVAL[NA.EST] <- NA
         NA
        PVAL <- NA
      }
      if (!all(NA.EST) & TTEST) {
        if (is.null(exact))
          exact <- (n < 50)
        EXACT.TIES <- exact & !TIES

        if (any(EXACT.TIES)) {
          exactq <- round((rhohat + 1) * n * (n - 1)/4)
          STATISTIC <- exactq
          pkendall <- function(q, n) .Call(C_pKendall, q, n)
          PVAL <- switch(alternative,
                         two.sided = {
                           loc.exactq <- exactq > n * (n - 1)/4
                           pv <- matrix(NA, nrow = p, ncol = q)
                           if (any(loc.exactq))
                             pv[loc.exactq] <- 1 - pkendall(exactq[loc.exactq] - 1, n)
                           if (any(!loc.exactq))
                             pv[!loc.exactq] <- pkendall(exactq[!loc.exactq], n)
                           min(2 * pv, 1)
                         },
                         greater = 1 - pkendall(exactq - 1, n),
                         less = pkendall(exactq, n))
        }



        ###???





        if (any(!EXACT.TIES)) {
          xties <- table(x[duplicated(x)]) + 1
          yties <- table(y[duplicated(y)]) + 1
          T0 <- n * (n - 1)/2
          T1 <- sum(xties * (xties - 1))/2
          T2 <- sum(yties * (yties - 1))/2
          S <- rhohat * sqrt((T0 - T1) * (T0 - T2))
          v0 <- n * (n - 1) * (2 * n + 5)
          vt <- sum(xties * (xties - 1) * (2 * xties +
                                             5))
          vu <- sum(yties * (yties - 1) * (2 * yties +
                                             5))
          v1 <- sum(xties * (xties - 1)) * sum(yties *
                                                 (yties - 1))
          v2 <- sum(xties * (xties - 1) * (xties - 2)) *
            sum(yties * (yties - 1) * (yties - 2))
          var_S <- (v0 - vt - vu)/18 + v1/(2 * n * (n -
                                                      1)) + v2/(9 * n * (n - 1) * (n - 2))
          if (exact && TIES)
            warning("Cannot compute exact p-value with ties")
          if (continuity)
            S <- sign(S) * (abs(S) - 1)
          STATISTIC <- c(z = S/sqrt(var_S))
          PVAL <- switch(alternative, less = pnorm(STATISTIC),
                         greater = pnorm(STATISTIC, lower.tail = FALSE),
                         two.sided = 2 * min(pnorm(STATISTIC), pnorm(STATISTIC,
                                                                     lower.tail = FALSE)))
        }
      }
    }
    else {
      method <- "Spearman's rank correlation rho"
      if (is.null(exact))
        exact <- TRUE
      names(NVAL) <- "rho"
      rhohat <- cor(rank(x), rank(y))
      ESTIMATE <- c(rho = rhohat)
      if (!is.finite(ESTIMATE)) {
        ESTIMATE[] <- NA
        STATISTIC <- c(S = NA)
        PVAL <- NA
      }
      else {
        pspearman <- function(q, n, lower.tail = TRUE) {
          if (n <= 1290 && exact)
            .Call(C_pRho, round(q) + 2 * lower.tail,
                  n, lower.tail)
          else {
            den <- (n * (n^2 - 1))/6
            if (continuity)
              den <- den + 1
            rhohat <- 1 - q/den
            pt(rhohat/sqrt((1 - rhohat^2)/(n - 2)), df = n - 2,
               lower.tail = !lower.tail)
          }
        }
        q <- (n^3 - n) * (1 - rhohat)/6
        STATISTIC <- c(S = q)
        if (TIES && exact) {
          exact <- FALSE
          warning("Cannot compute exact p-value with ties")
        }
        PVAL <- switch(alternative, two.sided = {
          p <- if (q > (n^3 - n)/6) pspearman(q, n, lower.tail = FALSE) else pspearman(q,
                                                                                       n, lower.tail = TRUE)
          min(2 * p, 1)
        }, greater = pspearman(q, n, lower.tail = TRUE),
        less = pspearman(q, n, lower.tail = FALSE))
      }
    }
  }


  RVAL <- list(statistic = STATISTIC, parameter = PARAMETER,
               p.value = as.numeric(PVAL), estimate = ESTIMATE, null.value = NVAL,
               alternative = alternative, method = method, data.name = DNAME)
  if (conf.int)
    RVAL <- c(RVAL, list(conf.int = cint))
  class(RVAL) <- "htest"
  RVAL
}


# A wrapper of 'qvalue::qvalue' to solve a bug in ???
# @importFrom qvalue qvalue
#
q.adjust <- function (p,
                      fdr.level = NULL,
                      pfdr = FALSE,
                      lfdr.out = TRUE,
                      pi0.method = c("smoother", "bootstrap")) {
  # lambda
  lambda <- seq(0.0, .9, .05)
  # call qvalue
  qvalues <- qvalue (p,
                     fdr.level = fdr.level,
                     pfdr = pfdr,
                     lfdr.out = lfdr.out,
                     lambda = lambda,
                     pi0.method = pi0.method)

  # Return q-values
  return(qvalues$qvalues)
}

## A self implementation of q-value correction for multiple comparison
## Storey JD. (2002) A direct approach to false discovery rates. Journal of the Royal Statistical Society, Series B, 64: 479-498.
## TO BE DEFINE
## TO BE WRITTEN
qvalue <- function (p, fdr.level = NULL, pfdr = FALSE,
                    lfdr.out = TRUE, lambda = NULL,
                    pi0 = NULL, pi0.method, ...) {
  return(p)
}
