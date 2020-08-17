
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
[![](https://img.shields.io/github/languages/code-size/villegar/MetaPipe.svg)](https://github.com/villegar/MetaPipe)
[![](https://www.r-pkg.org/badges/version/MetaPipe?color=red)](https://cran.r-project.org/package=MetaPipe)
[![](https://codecov.io/gh/villegar/MetaPipe/branch/master/graph/badge.svg)](https://codecov.io/gh/villegar/MetaPipe)
<!-- [![Dependencies](https://tinyverse.netlify.com/badge/MetaPipe)](https://cran.r-project.org/package=MetaPipe) -->
<!-- [![CRAN checks](https://cranchecks.info/badges/summary/MetaPipe)](https://cran.r-project.org/web/checks/check_results_MetaPipe.html) -->

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
# install.packages("remotes")
remotes::install_github("villegar/MetaPipe")
```

## Example

<!-- This is a basic example which shows you how to solve a common problem: -->

You should start by loading `MetaPipe` on your session.

``` r
library(MetaPipe)
```

### Loading raw data

`MetaPipe` only accepts comma-separated values (CSV) files with the
following
structure:

<table class="table table-striped table-hover" style="width: auto !important; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:center;color: #EEEEEE !important;background-color: #363B74 !important;vertical-align: middle;">

ID

</th>

<th style="text-align:center;color: #EEEEEE !important;background-color: #363B74 !important;vertical-align: middle;">

\[Property\]<sub>1</sub>

</th>

<th style="text-align:center;color: #EEEEEE !important;background-color: #363B74 !important;vertical-align: middle;">

…

</th>

<th style="text-align:center;color: #EEEEEE !important;background-color: #363B74 !important;vertical-align: middle;">

\[Property\]<sub>M</sub>

</th>

<th style="text-align:center;color: #EEEEEE !important;background-color: #363B74 !important;vertical-align: middle;">

\[Trait\]<sub>1</sub>

</th>

<th style="text-align:center;color: #EEEEEE !important;background-color: #363B74 !important;vertical-align: middle;">

…

</th>

<th style="text-align:center;color: #EEEEEE !important;background-color: #363B74 !important;vertical-align: middle;">

\[Trait\]<sub>N</sub>

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:center;">

 

</td>

<td style="text-align:center;">

</td>

<td style="text-align:center;">

</td>

<td style="text-align:center;">

</td>

<td style="text-align:center;">

</td>

<td style="text-align:center;">

</td>

<td style="text-align:center;">

</td>

</tr>

<tr>

<td style="text-align:center;">

 

</td>

<td style="text-align:center;">

</td>

<td style="text-align:center;">

</td>

<td style="text-align:center;">

</td>

<td style="text-align:center;">

</td>

<td style="text-align:center;">

</td>

<td style="text-align:center;">

</td>

</tr>

</tbody>

</table>

Where the first column (`ID`) should be an unique identifier for each
entry, if there are repeated values `MetaPipe` will aggregate (mean) and
replace them by a single row. The data structure can have 0 to M
properties, including categorical and numeric. Finally, at least one one
trait is expected.

The function call is as follows:

``` r
load_raw(raw_data_filename = "FILE.CSV", excluded_columns = c(...))
```

where `raw_data_filename` is the filename containg the raw data, either
absolute or relative paths are accepted. Next, the argument
`excluded_columns` is a vector containing the indices of the properties,
e.g. `c(2, 3, ..., M)`.

``` r
# Dummy data set
example_data <- data.frame(ID = c(1,2,3,4,5),
                           P1 = c("one", "two", "three", "four", "five"), 
                           F1 = rnorm(5), 
                           F2 = rnorm(5))
## Write to disk
write.csv(example_data, "example_data.csv", row.names = FALSE)
## Create copy with duplicated rows
write.csv(example_data[c(1:5, 1, 2), ], "example_data_dup.csv", row.names = FALSE)

# Load the data
load_raw("example_data.csv", c(2))
#>   ID    P1         F1         F2
#> 1  1   one  0.7385226 -0.5075755
#> 2  2   two  1.9370058  0.9575687
#> 3  3 three -0.6704865 -0.3290465
#> 4  4  four -0.4076667  1.2337800
#> 5  5  five -1.0915699  0.5501291
load_raw("example_data_dup.csv", c(2))
#>   ID    P1         F1         F2
#> 1  1   one  0.7385226 -0.5075755
#> 2  2   two  1.9370058  0.9575687
#> 3  3 three -0.6704865 -0.3290465
#> 4  4  four -0.4076667  1.2337800
#> 5  5  five -1.0915699  0.5501291
```

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
