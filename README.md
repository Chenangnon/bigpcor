
# bigpcor
The package `bigpcor` is a `R` library for a fast computation of partial correlations conditioning on a large set of confounders.

The library takes advantage of linear algebra to reduce the number of actual covariance matrix inversions required to compute pairwise partial correlations given a large set of variables to condition on.

## Installation

You can install the development version of bigpcor from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Chenangnon/bigpcor")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(bigpcor)
## basic example code
```
