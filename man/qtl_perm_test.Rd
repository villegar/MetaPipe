% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metapipe.R
\name{qtl_perm_test}
\alias{qtl_perm_test}
\title{QTL mapping permutation test}
\usage{
qtl_perm_test(
  x_data,
  cpus = 1,
  qtl_method = "par-scanone",
  raw_data_normalised = NULL,
  lod_threshold = 3,
  parametric = TRUE,
  n_perm = 1000,
  plots_dir = tempdir(),
  ...
)
}
\arguments{
\item{x_data}{Cross-data frame containing genetic map data and traits.}

\item{cpus}{Number of CPUs to be used in the computation.}

\item{qtl_method}{QTL mapping method.}

\item{raw_data_normalised}{Normalised raw data, see 
\code{\link{assess_normality}}.}

\item{lod_threshold}{LOD score threshold to look up for significant QTLs}

\item{parametric}{Boolean flag to indicate whether or not \code{x_data} 
contains parametric (normal) traits.}

\item{n_perm}{Number of permutations.}

\item{plots_dir}{Output directory for plots.}

\item{...}{Arguments passed on to 
\code{\link[qtl:scanone]{qtl::scanone}}.}
}
\value{
Data frame containing the significant QTLs information.
}
\description{
Perform a QTL mapping permutation test using the 
\code{\link[qtl:scanone]{qtl:scanone(...)}} function to find significant QTL.
}
\examples{
\donttest{
# Create temp dir
tmp <- tempdir()
dir.create(tmp, showWarnings = FALSE, recursive = TRUE)

# Toy dataset
excluded_columns <- c(1, 2)
population <- 5
seed <- 123
set.seed(seed)
example_data <- data.frame(ID = 1:population,
                           P1 = c("one", "two", "three", "four", "five"),
                           T1 = rnorm(population),
                           T2 = rnorm(population))

output <- MetaPipe::assess_normality(example_data, 
                                     excluded_columns, 
                                     show_stats = FALSE,
                                     out_prefix = paste0(tmp, "/tmp"))

# Create and store random genetic map (for testing only)
genetic_map <- MetaPipe:::random_map(population = population, 
                                     seed = seed)

# Load cross file with genetic map and raw data for normal traits
x <- MetaPipe::read.cross(genetic_map, output$norm)

x <- qtl::calc.genoprob(x, step = 1, error.prob = 0.001)
x_scone <- MetaPipe::qtl_scone(x, 1, model = "normal", method = "hk")
x_qtl_perm <- MetaPipe::qtl_perm_test(x, 
                                      n_perm = 5, 
                                      model = "normal", 
                                      method = "hk",
                                      plots_dir = tmp)
x_qtl_perm_1000 <- MetaPipe::qtl_perm_test(x, 
                                           n_perm = 1000, 
                                           model = "normal", 
                                           method = "hk",
                                           plots_dir = tmp)

# F1 Seedling Ionomics dataset
data(ionomics) # Includes some missing data
data(father_riparia) # Genetic map
ionomics_rev <- MetaPipe::replace_missing(ionomics, 
                                          excluded_columns = c(1, 2),
                                          replace_na =  TRUE,
                                          out_prefix = paste0(tmp, "/tmp"))
ionomics_normalised <- 
  MetaPipe::assess_normality(ionomics_rev,
                             excluded_columns = c(1, 2),
                             out_prefix = file.path(tmp, "ionomics"),
                             transf_vals = c(2, exp(1)),
                             show_stats = FALSE)

# Load cross file with genetic map and raw data for normal traits
x <- MetaPipe::read.cross(father_riparia, 
                          ionomics_normalised$norm,
                          genotypes = c("nn", "np", "--"))
                          
set.seed(seed)
x <- qtl::jittermap(x)
x <- qtl::calc.genoprob(x, step = 1, error.prob = 0.001)

x_scone <- MetaPipe::qtl_scone(x, 1, model = "normal", method = "hk")
x_qtl_perm <- MetaPipe::qtl_perm_test(x, 
                                      n_perm = 5, 
                                      model = "normal", 
                                      method = "hk",
                                      plots_dir = tmp)
x_qtl_perm_1000 <- MetaPipe::qtl_perm_test(x, 
                                           n_perm = 1000, 
                                           model = "normal", 
                                           method = "hk",
                                           plots_dir = tmp)

# Clean temporal directory
# unlink(tmp, recursive = TRUE, force = TRUE)
MetaPipe:::tidy_up(tmp)
}

}
\seealso{
Other QTL mapping functions: 
\code{\link{qtl_scone}()},
\code{\link{read.cross}()}
}
\concept{QTL mapping functions}
