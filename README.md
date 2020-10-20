
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- Extra CSS -->

<!-- Utilitary functions -->

# MetaPipe <img src="https://raw.githubusercontent.com/villegar/MetaPipe/master/inst/images/metapipe.png" alt="logo" align="right" height=200px/>

MetaPipe: A High-Performance Computing Pipeline for QTL Mapping of Large
Ionomic and Metabolomic Datasets <!-- badges: start -->
<!-- [![Build Status](https://travis-ci.com/villegar/MetaPipe.svg?branch=master)](https://travis-ci.com/villegar/MetaPipe) -->
<!-- [![](https://travis-ci.com/villegar/MetaPipe.svg?branch=master)](https://travis-ci.com/villegar/MetaPipe) -->
[![](https://www.r-pkg.org/badges/version/MetaPipe?color=red)](https://cran.r-project.org/package=MetaPipe)
[![](https://img.shields.io/badge/devel%20version-0.0.0.900-yellow.svg)](https://github.com/villegar/MetaPipe)
[![](https://codecov.io/gh/villegar/MetaPipe/branch/master/graph/badge.svg)](https://codecov.io/gh/villegar/MetaPipe)
[![R build
status](https://github.com/villegar/MetaPipe/workflows/R-CMD-check/badge.svg)](https://github.com/villegar/MetaPipe/actions)
<!-- [![](https://img.shields.io/github/languages/code-size/villegar/MetaPipe.svg)](https://github.com/villegar/MetaPipe) -->
<!-- [![Dependencies](https://tinyverse.netlify.com/badge/MetaPipe)](https://cran.r-project.org/package=MetaPipe) -->
<!-- [![CRAN checks](https://cranchecks.info/badges/summary/MetaPipe)](https://cran.r-project.org/web/checks/check_results_MetaPipe.html) -->

<!-- badges: end -->

## Overview

The goal of MetaPipe is to provide an easy to use and powerful tool
capable of performing QTL mapping analyses.
<!-- on metabolomics data. -->

## Installation

You can install the released version of MetaPipe from
[CRAN](https://cran.r-project.org/package=MetaPipe) with:

``` r
install.packages("MetaPipe")
```

And the development version from
[GitHub](https://github.com/villegar/MetaPipe) with:

``` r
# install.packages(c("hexSticker", "kableExtra", "qpdf", "remotes")
remotes::install_github("villegar/MetaPipe", build_vignettes = TRUE)
```

## Example

<!-- This is a basic example which shows you how to solve a common problem: -->

<!-- You should start by loading `MetaPipe` on your session. -->

<!-- ```{r example} -->

<!-- library(MetaPipe) -->

<!-- ``` -->

### Load raw data

For details about the data structure and extended documentation, see the
vignette [Load Raw
Data](https://villegar.github.io/MetaPipe/articles/load-raw-data).

``` r
vignette("load-raw-data", package = "MetaPipe")
```

#### Function call

``` r
MetaPipe::load_raw(raw_data_filename = "FILE.CSV", excluded_columns = c(...))
```

where `raw_data_filename` is the filename containing the raw data, both
absolute and relative paths are accepted. Next, the argument
`excluded_columns` is a vector containing the indices of the properties,
e.g. `c(2, 3, ..., M)`.

``` r
# F1 Seedling Ionomics dataset
ionomics_path <- system.file("extdata", 
                             "ionomics.csv", 
                             package = "MetaPipe", 
                             mustWork = TRUE)
ionomics <- MetaPipe::load_raw(ionomics_path)
knitr::kable(ionomics[1:5, 1:8])
```

| ID     | SampleWeight |     Ca44 |      K39 |      P31 |       Li7 |      B11 |      Na23 |
| :----- | -----------: | -------: | -------: | -------: | --------: | -------: | --------: |
| E\_199 |           79 | 32675.79 | 6051.023 | 2679.338 | 0.1159068 | 23.32975 |  9.372606 |
| E\_209 |           81 | 28467.95 | 5642.651 | 2075.403 | 0.0104801 | 27.31206 |  8.787553 |
| E\_035 |           81 | 27901.35 | 7357.856 | 2632.343 | 0.0561879 | 16.87480 | 14.369062 |
| E\_197 |           79 | 27855.36 | 5225.275 | 1761.725 | 0.0104453 | 25.34740 | 11.009597 |
| E\_016 |           79 | 27377.40 | 6141.001 | 2145.715 | 0.0172996 | 24.64500 |  6.999958 |

### Replace missing data

For extended documentation, see the vignette [Replace Missing
Data](https://villegar.github.io/MetaPipe/articles/replace-missing-data).

``` r
vignette("replace-missing-data", package = "MetaPipe")
```

#### Function call

``` r
MetaPipe::replace_missing(raw_data = example_data, 
                          excluded_columns = c(2), 
                          # Optional
                          out_prefix = "metapipe", 
                          prop_na = 0.5, 
                          replace_na = FALSE)
```

where `raw_data` is a data frame containing the raw data, as described
in [Load Raw Data](#load-raw-data) and `excluded_columns` is a vector
containing the indices of the properties, e.g. `c(2, 3, ..., M)`. The
other arguments are optional, `out_prefix` is the prefix for output
files, `prop_na` is the proportion of NA values (used to drop traits),
and `replace_na` is a logical flag to indicate whether or not `NA`s
should be replace by half of the minimum value within each variable.

``` r
# F1 Seedling Ionomics dataset
data(ionomics) # Includes some missing data
ionomics_rev <- MetaPipe::replace_missing(ionomics, c(1, 2))
ionomics_rev <- MetaPipe::replace_missing(ionomics, 
                                          excluded_columns = c(1, 2), 
                                          prop_na =  0.025)
#> The following trait was dropped because it has 2.5% or more missing values: 
#>  - Se78
ionomics_rev <- MetaPipe::replace_missing(ionomics, 
                                          excluded_columns = c(1, 2),
                                          replace_na =  TRUE)
knitr::kable(ionomics_rev[1:5, 1:8])
```

| ID     | SampleWeight |     Ca44 |      K39 |      P31 |       Li7 |      B11 |      Na23 |
| :----- | -----------: | -------: | -------: | -------: | --------: | -------: | --------: |
| E\_001 |           79 | 15894.22 | 5888.311 | 1743.118 | 0.0128699 | 18.66673 |  6.970224 |
| E\_002 |           93 | 13155.45 | 7013.400 | 2244.684 | 0.0119316 | 14.47693 |  5.866392 |
| E\_004 |           97 | 14182.51 | 7966.273 | 2311.057 | 0.0212316 | 14.71313 | 10.251955 |
| E\_005 |           82 | 22550.82 | 7514.089 | 2315.675 | 0.0233063 | 20.10630 | 11.773697 |
| E\_006 |           99 | 15982.76 | 7608.464 | 1995.193 | 0.0588128 | 12.97801 | 11.043837 |

### Assess normality

For extended documentation, see the vignette [Assess
Normality](https://villegar.github.io/MetaPipe/articles/assess-normality).

``` r
vignette("assess-normality", package = "MetaPipe")
```

`MetaPipe` assesses the normality of variables (traits) by performing a
Shapiro-Wilk test on the raw data (see [Load Raw
Data](https://villegar.github.io/MetaPipe/articles/load-raw-data.html)
and [Replace Missing
Data](https://villegar.github.io/MetaPipe/articles/replace-missing-data.html).
Based on whether or not the data approximates a normal distribution, an
array of transformations will be computed, and the normality assessed
one more time.

#### Function call

``` r
MetaPipe::assess_normality(raw_data = raw_data, 
                           excluded_columns = c(2, 3, ..., M), 
                           # Optional
                           cpus = 1, 
                           out_prefix = "metapipe", 
                           plots_dir = getwd(), 
                           transf_vals = c(2, exp(1), 3, 4, 5, 6, 7, 8, 9, 10),
                           alpha = 0.05,
                           pareto_scaling = FALSE,
                           show_stats = TRUE)
```

where `raw_data` is a data frame containing the raw data, as described
in [Load Raw
Data](https://villegar.github.io/MetaPipe/articles/load-raw-data.html)
and `excluded_columns` is a vector containing the indices of the
properties, e.g. `c(2, 3, ..., M)`. The other arguments are optional,
`cpus` is the number of cores to use, in other words, the number of
concurrent traits to process, `out_prefix` is the prefix for output
files, `plots_dir` is the output directory where the plots will be
stored, `transf_vals` is a vector containing the transformation values
to be used when transforming the original data, `alpha` is the
significance level for the *Wilk-Shapiro* tests, `pareto_scaling` is a
boolean flag to indicate whether or not to scale the traits to the same
scale, and `show_stats` is a boolean flag to show or hide some general
statistics of the normalisation process.

``` r
# F1 Seedling Ionomics dataset
data(ionomics) # Includes some missing data
ionomics_rev <- MetaPipe::replace_missing(ionomics, 
                                          excluded_columns = c(1, 2),
                                          replace_na =  TRUE)
ionomics_normalised <- 
  MetaPipe::assess_normality(ionomics_rev,
                             excluded_columns = c(1, 2),
                             transf_vals = c(2, exp(1)),
                             out_prefix = "README-ionomics",
                             plots_dir = "man/figures/",
                             pareto_scaling = FALSE)
#> Total traits (excluding all NAs traits):     21
#> Normal traits (without transformation):      2
#> Normal traits (transformed):                 4
#> Total normal traits:                         6
#> Total skewed traits:                         15
#> 
#> Transformations summary:
#>  f(x)      Value     # traits  
#>  log       2         3         
#>  root      e         1

# Extract normalised features
ionomics_norm <- ionomics_normalised$norm
ionomics_skew <- ionomics_normalised$skew
```

The function call to `MetaPipe::assess_normality` will print a summary
of the transformations performed (if any), as well as an overview of the
number of traits that should be considered *normal* and *skewed*. Next,
we can preview some of the partial output of the normality assessment
process:

``` r
# Normal traits
knitr::kable(ionomics_norm[1:5, ]) 
```

| ID     |     Ca44 |      B11 |     Na23 |     Mg26 |     Rb85 |      Sr88 |
| :----- | -------: | -------: | -------: | -------: | -------: | --------: |
| E\_001 | 15894.22 | 4.222397 | 2.042740 | 10.77021 | 1.555742 |  7.347059 |
| E\_002 | 13155.45 | 3.855684 | 1.917202 | 10.54095 | 2.058711 |  6.890243 |
| E\_004 | 14182.51 | 3.879033 | 2.354263 | 10.51931 | 2.198422 |  9.025915 |
| E\_005 | 22550.82 | 4.329576 | 2.477233 | 11.13450 | 1.791578 | 15.292360 |
| E\_006 | 15982.76 | 3.697997 | 2.419593 | 11.72734 | 2.229866 | 13.901449 |

``` r
# Skewed traits (partial output)
knitr::kable(ionomics_skew[1:5, 1:8])
```

| ID     |      K39 |      P31 |       Li7 |      Al27 |      S34 |     Fe54 |     Mn55 |
| :----- | -------: | -------: | --------: | --------: | -------: | -------: | -------: |
| E\_001 | 5888.311 | 1743.118 | 0.0128699 |  3.845879 | 1152.944 | 27.59340 | 54.53991 |
| E\_002 | 7013.400 | 2244.684 | 0.0119316 |  5.825639 | 1600.442 | 35.49159 | 52.57114 |
| E\_004 | 7966.273 | 2311.057 | 0.0212316 |  8.036047 | 1039.098 | 39.13434 | 36.66475 |
| E\_005 | 7514.089 | 2315.675 | 0.0233063 |  9.482051 | 1091.607 | 40.22041 | 43.24368 |
| E\_006 | 7608.464 | 1995.193 | 0.0588128 | 29.329605 | 1096.871 | 75.23614 | 53.64705 |

Among the transformed traits, we have `B11` and `Na23`. Both of which
seem to be skewed, but after a simple transformation, can be classify as
normalised traits.

<img src="man/figures/HIST_5_LOG_2_B11.png" width="45%" />
<img src="man/figures/HIST_6_ROOT_e_Na23.png" width="45%" />

### QTL mapping
