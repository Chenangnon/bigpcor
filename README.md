
## bigpcor
The package `bigpcor` is a `R` library for a fast computation of partial correlations conditioning on a large set of confounders.

The library takes advantage of linear algebra to reduce the number of actual covariance matrix inversions required to compute pairwise partial correlations given a large set of variables to condition on.

The library depends on R packages: `Rdpack`, `methods`, `parallel`, and `propagate`.

## Installation

You can install the development version of bigpcor from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Chenangnon/bigpcor")
```

## Example

This is a basic example which shows you how to solve a common problem.
The example mimicks the selection of confounders of gene expressions in a genomic network. We consider a simulated data with 50 genes ('fiftygenes').

#* Task: compute pairwise partial correlations between expressions on one hand and candidate confouders on the other hand, conditioning on genetic variants.

First load the package
``` r
library(bigpcor)
```

#* Set number of genetic variants (q), genes (p) and confounders (u) in the dataset 'fiftygenes'

``` r
p <- q <- 50
u <- 100
```

#* Using 'bigpcor'
```r
Time <- system.time({
  res <- bigpcor (x = fiftygenes[, (q+1):(q+p)],
                  y = fiftygenes[, (q+p+1):(q+p+u)],
                  z = fiftygenes[, 1:q])
})
```

A glance at the results for the first five genes and five confounders
```r
res[1:5, 1:5]

# Elapsed time (about 2 seconds)
Time
```

Library 'ppcor' is required to run the remaining of the example


#* Comparing 'bigpcor' with a repeated call to 'pcor' of library "ppcor"
# Install 'ppcor' if not installed
```r
if (!"ppcor" %in% rownames(installed.packages()))
  install.packages('ppcor', dependencies = TRUE)
require(ppcor)
```

Define a function to disable the display of warnings from 'ppcor::pcor'. This is Martin Maechler's 'tryCatch.W.E' which catches and saves both errors and warnings, and in the case of a warning, also keeps the computed result.
```r
catch.conditions <- function (expr) {
  W <- NULL
  w.handler <- function(w) {
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                  warning = w.handler), warning = W)
}
```

#* Use the 'pcor' function of library "ppcor"
The following took about 415 seconds (~7 minutes) on a MAC OS Ventura 13.4.1 system with 16 GB RAM
```r
W.E. <- catch.conditions({
  Time0 <- system.time({
    res0 <- apply(expand.grid(E = (q+1):(q+p), C = (q+p+1):(q+p+u)),# Combinations of E and C
                  MARGIN = 1,
                  FUN = function(x) { # call 'ppcor::pcor'
                    pcor (fiftygenes[, c(x, 1:q)])$estimate[1,2]
                  })
  })
})
```

Run the following to display catched warning(s)
W.E.$warning

Organize res0 into a matrix (rows for genes and columns for confounders)
```r
res0 <- matrix(res0, nrow = p, ncol = u, byrow = FALSE)
```

A glance at the results for the first five genes and five confounders
```r
res0[1:5, 1:5]
```

Compare results to 10 decimal places
```r
sum(round(res, 10) != round(res0, 10)) # zero means no difference
```

```r
Time0
Time0[3]/Time[3] # 'bigpcor' is about 190 fold faster than 'ppcor' in this example
```
