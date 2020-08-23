
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- Extra CSS -->

<!-- Utilitary functions -->

<!-- /Library/Frameworks/R.framework/Versions/4.0/Resources/library/MetaPipe/images/metapipe.png -->

# MetaPipe <img src="https://raw.githubusercontent.com/villegar/MetaPipe/master/inst/images/metapipe.png" alt="metapipe-logo" align="right" height=200px/>

MetaPipe: A High-Performance Computing pipeline for QTL mapping of large
metabolomic datasets <!-- badges: start -->
<!-- [![Build Status](https://travis-ci.com/villegar/MetaPipe.svg?branch=master)](https://travis-ci.com/villegar/MetaPipe) -->
[![](https://travis-ci.com/villegar/MetaPipe.svg?branch=master)](https://travis-ci.com/villegar/MetaPipe)
[![](https://img.shields.io/badge/devel%20version-0.0.1-blue.svg)](https://github.com/villegar/MetaPipe)
[![](https://img.shields.io/github/languages/code-size/villegar/MetaPipe.svg)](https://github.com/villegar/MetaPipe)<!-- [![](https://www.r-pkg.org/badges/version/MetaPipe?color=red)](https://cran.r-project.org/package=MetaPipe) -->
[![](https://codecov.io/gh/villegar/MetaPipe/branch/master/graph/badge.svg)](https://codecov.io/gh/villegar/MetaPipe)
[![R build
status](https://github.com/villegar/MetaPipe/workflows/R-CMD-check/badge.svg)](https://github.com/villegar/MetaPipe/actions)
<!-- [![Dependencies](https://tinyverse.netlify.com/badge/MetaPipe)](https://cran.r-project.org/package=MetaPipe) -->
<!-- [![CRAN checks](https://cranchecks.info/badges/summary/MetaPipe)](https://cran.r-project.org/web/checks/check_results_MetaPipe.html) -->
<!-- [![R build status](https://github.com/villegar/MetaPipe/workflows/R-CMD-check/badge.svg)](https://github.com/villegar/MetaPipe/actions) -->

<!-- badges: end -->

## Overview

The goal of MetaPipe is to provide an easy to use and powerful tool
capable of performing QTL mapping analyses.
<!-- on metabolomics data. -->

## Installation

<!-- You can install the released version of MetaPipe from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("MetaPipe") -->

<!-- ``` -->

<!-- And the development version from [GitHub](https://github.com/) with: -->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages(c("hexSticker", "kableExtra", "qpdf", remotes")
remotes::install_github("villegar/MetaPipe", build_vignettes = TRUE)
```

## Example

<!-- This is a basic example which shows you how to solve a common problem: -->

You should start by loading `MetaPipe` on your session.

``` r
library(MetaPipe)
```

### Load raw data

For details about the data structure and extended documentation, see the
vignette [Load Raw
Data](https://villegar.github.io/MetaPipe/articles/load-raw-data).

``` r
vignette("load-raw-data", package = "MetaPipe")
```

#### Function call

``` r
load_raw(raw_data_filename = "FILE.CSV", excluded_columns = c(...))
```

where `raw_data_filename` is the filename containg the raw data, both
absolute and relative paths are accepted. Next, the argument
`excluded_columns` is a vector containing the indices of the properties,
e.g. `c(2, 3, ..., M)`.

``` r
# Toy dataset
set.seed(123)
example_data <- data.frame(ID = c(1,2,3,4,5),
                           P1 = c("one", "two", "three", "four", "five"), 
                           T1 = rnorm(5), 
                           T2 = rnorm(5),
                           T3 = c(NA, rnorm(4)),                     #  20 % NAs
                           T4 = c(NA, 1.2, -0.5, NA, 0.87),          #  40 % NAs
                           T5 = NA)                                  # 100 % NAs
## Write to disk
write.csv(example_data, "example_data.csv", row.names = FALSE)

# Load the data
load_raw("example_data.csv", c(2))
#>   ID    P1          T1         T2        T3    T4 T5
#> 1  1   one -0.56047565  1.7150650        NA    NA NA
#> 2  2   two -0.23017749  0.4609162 1.2240818  1.20 NA
#> 3  3 three  1.55870831 -1.2650612 0.3598138 -0.50 NA
#> 4  4  four  0.07050839 -0.6868529 0.4007715    NA NA
#> 5  5  five  0.12928774 -0.4456620 0.1106827  0.87 NA
```

### Replace missing data

For extended documentation, see the vignette [Replace Missing
Data](https://villegar.github.io/MetaPipe/articles/replace-missing-data).

``` r
vignette("replace-missing-data", package = "MetaPipe")
```

#### Function call

``` r
replace_missing(raw_data = example_data, 
                excluded_columns = c(2), 
                # Optional
                out_prefix = "metapipe", prop_na = 0.5, replace_na = FALSE)
```

where `raw_data` is a data frame containing the raw data, as described
in [Load Raw Data](#load-raw-data) and `excluded_columns` is a vector
containing the indices of the properties, e.g. `c(2, 3, ..., M)`. The
other arguments are optional, `out_prefix` is the prefix for output
files, `prop_na` is the proportion of NA values (used to drop traits),
and `replace_na` is a logical flag to indicate whether or not NAs should
be replace by half of the minimum value.

``` r
# Inserting missing values manually
example_data$T1[2:3] <- NA
example_data$T2[4] <- NA

# No changes expected
replace_missing(example_data, c(2))
#> The following traits were dropped because they have 50% or more missing values: 
#> T5
#>   ID    P1          T1         T2        T3    T4
#> 1  1   one -0.56047565  1.7150650        NA    NA
#> 2  2   two          NA  0.4609162 1.2240818  1.20
#> 3  3 three          NA -1.2650612 0.3598138 -0.50
#> 4  4  four  0.07050839         NA 0.4007715    NA
#> 5  5  five  0.12928774 -0.4456620 0.1106827  0.87

# Traits with 25% of NA should be dropped
replace_missing(example_data, c(2), prop_na =  0.25)
#> The following traits were dropped because they have 25% or more missing values: 
#> T1, T4, T5
#>   ID    P1         T2        T3
#> 1  1   one  1.7150650        NA
#> 2  2   two  0.4609162 1.2240818
#> 3  3 three -1.2650612 0.3598138
#> 4  4  four         NA 0.4007715
#> 5  5  five -0.4456620 0.1106827

# NAs should be replaced by half of the minimum value
replace_missing(example_data, c(2), replace_na =  TRUE)
#> The following traits were dropped because they have 100% missing values: 
#> T5
#>   ID    P1          T1         T2         T3    T4
#> 1  1   one -0.56047565  1.7150650 0.05534136 -0.25
#> 2  2   two -0.28023782  0.4609162 1.22408180  1.20
#> 3  3 three -0.28023782 -1.2650612 0.35981383 -0.50
#> 4  4  four  0.07050839 -0.6325306 0.40077145 -0.25
#> 5  5  five  0.12928774 -0.4456620 0.11068272  0.87
```

### Assess normality

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->

<!-- ```{r cars} -->

<!-- summary(cars) -->

<!-- ``` -->

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. -->

<!-- You can also embed plots, for example: -->

<!-- ```{r pressure, echo = FALSE} -->

<!-- plot(pressure) -->

<!-- ``` -->

<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub! -->
