#' @keywords internal
"_PACKAGE"

#' Load raw data
#' 
#' Load raw data from disk and aggregate (using the \code{mean} function) 
#' observations with duplicated IDs (first column). Non-numeric columns must
#' be excluded using the \code{excluded_columns} parameter.
#' 
#' @importFrom stats aggregate
#' @importFrom stats na.omit
#' @importFrom utils read.csv
#' @importFrom utils write.csv
#' 
#' @param raw_data_filename Filename containing the raw data, it can be a
#'     relative path (e.g. \code{"my_input.csv"}) or an absolute path (e.g. 
#'     \code{"/path/to/my_input.csv"}).
#' @param excluded_columns Numeric vector containing the indices of the dataset 
#'     properties that are non-numeric, excluded columns.
#'
#' @return Data frame with the pre-processed raw data.
#' @export
#'
#' @examples
#' # Toy dataset
#' example_data <- data.frame(ID = c(1,2,3,4,5), 
#'                            P1 = c("one", "two", "three", "four", "five"), 
#'                            T1 = rnorm(5), 
#'                            T2 = rnorm(5))
#' write.csv(example_data, "example_data.csv", row.names = FALSE)
#' write.csv(example_data[c(1:5, 1, 2), ], 
#'           "example_data_dup.csv", 
#'           row.names = FALSE)
#' knitr::kable(MetaPipe::load_raw("example_data.csv", c(1, 2)))
#' knitr::kable(MetaPipe::load_raw("example_data_dup.csv", c(1, 2)))
#' 
#' 
#' # F1 Seedling Ionomics dataset
#' ionomics_path <- system.file("extdata", 
#'                              "ionomics.csv", 
#'                              package = "MetaPipe", 
#'                              mustWork = TRUE)
#' ionomics <- MetaPipe::load_raw(ionomics_path)
#' knitr::kable(ionomics[1:5, 1:8])
#' 
#' # Clean up example outputs
#' MetaPipe:::tidy_up("example_data")
load_raw <- function(raw_data_filename, excluded_columns = NULL) {
  # Load and clean raw data
  raw_data <- read.csv(raw_data_filename, stringsAsFactors = FALSE)
  
  # Exclude first column, ID
  excluded_columns <- check_types(raw_data, unique(c(1, excluded_columns)))
  
  # Enforce first column name, ID
  colnames(raw_data)[1] <- "ID"
  
  # Aggregate data by row, computing the mean
  mean_raw_data <- aggregate(raw_data[, -excluded_columns],
                             by = list(raw_data[, 1]),
                             mean, 
                             na.action = na.omit)
  colnames(mean_raw_data)[1] <- colnames(raw_data)[1]
  
  # Joins the agreggated data with the excluded columns (properties)
  mean_raw_data <- dplyr::left_join(raw_data[, excluded_columns, drop = FALSE],
                                    mean_raw_data,
                                    by = colnames(raw_data)[1])
  
  # Find duplicated row values on the first column (ID)
  mean_raw_data <- mean_raw_data[!duplicated(mean_raw_data[, 1]), ]
  rownames(mean_raw_data) <- 1:nrow(mean_raw_data)
  return(mean_raw_data)
}

#' Replace missing values (\code{NA}s)
#' 
#' @description 
#' Replace missing values (\code{NA}s) in a dataset, the user can choose 
#' between two actions to handle missing data:
#' \enumerate{
#'   \item Drop traits (variables) that exceed a given threshold, 
#'   \code{prop_na}, a rate of missing (\code{NA}) and total observations.
#' 
#'   \item Replace missing values by half of the minimum within each trait.
#' }
#' 
#' Finally, if there are traits for which all entries are missing, these will 
#' be removed from the dataset and stored in a external CSV file called
#' \code{"<out_prefix>_NA_raw_data.csv"}.
#' 
#' @param raw_data Data frame containing the raw data.
#' @param excluded_columns Numeric vector containing the indices of the dataset 
#'     properties that are non-numeric, excluded columns.
#' @param out_prefix Prefix for output files and plots.
#' @param prop_na Proportion of missing/total observations, if a trait exceeds 
#'     this threshold and \code{replace_na = FALSE}, then it will be 
#'     dropped out.
#' @param replace_na Boolean flag to indicate whether or not missing values 
#'     should be replaced by half of the minimum value within each trait.
#'
#' @return Data frame containing the raw data without missing values.
#' @export
#'
#' @examples                                         
#' # Toy dataset                                        
#' example_data <- data.frame(ID = c(1,2,3,4,5), 
#'                            P1 = c("one", "two", "three", "four", "five"), 
#'                            T1 = rnorm(5), 
#'                            T2 = rnorm(5),
#'                            T3 = c(NA, rnorm(4)),                  #  20 % NAs
#'                            T4 = c(NA, 1.2, -0.5, NA, 0.87),       #  40 % NAs
#'                            T5 = NA)                               # 100 % NAs
#' MetaPipe::replace_missing(example_data, c(1, 2))
#' MetaPipe::replace_missing(example_data, c(1, 2), prop_na =  0.25)
#' MetaPipe::replace_missing(example_data, c(1, 2), replace_na =  TRUE)
#' 
#' 
#' # F1 Seedling Ionomics dataset
#' data(ionomics) # Includes some missing data
#' ionomics_rev <- MetaPipe::replace_missing(ionomics, c(1, 2))
#' ionomics_rev <- MetaPipe::replace_missing(ionomics, 
#'                                           excluded_columns = c(1, 2), 
#'                                           prop_na =  0.025)
#' ionomics_rev <- MetaPipe::replace_missing(ionomics, 
#'                                           excluded_columns = c(1, 2),
#'                                           replace_na =  TRUE)
#' knitr::kable(ionomics_rev[1:5, 1:8])
#' 
#' # Clean up example outputs
#' MetaPipe:::tidy_up("metapipe_")
replace_missing <- function(raw_data,
                            excluded_columns = NULL,
                            out_prefix = "metapipe",
                            prop_na = 0.5,
                            replace_na = FALSE) {
  # Exclude column 1, ID
  excluded_columns <- unique(c(1, excluded_columns))
  # excluded_columns <- check_types(raw_data, 
  #                                 unique(c(1, excluded_columns)))
  
  # Replace missing values by half of the minimum non-zero value for each trait.
  if (replace_na) {
    raw_data[, -excluded_columns] <- sapply(raw_data[, -excluded_columns], 
                                            rplc_na)
  } else {
    # Find which variables exceed the proportion of NAs threshold, prop_na
    idx <- which(colMeans(is.na(raw_data[, -excluded_columns])) >= prop_na)
    # Extract the indices from the original data
    idx <- colnames(raw_data) %in% names(idx)
    if(sum(idx) > 0) {
      # Store the dropped traits in a CSV file
      write.csv(raw_data[, c(excluded_columns, idx)],
                file = paste0(out_prefix, "_NA_raw_data.csv"),
                row.names = FALSE)
      msg <- paste0("The following trait",
                    ifelse(sum(idx) > 1, "s were ", " was "),
                    "dropped because ",
                    ifelse(sum(idx) > 1, "they have ", "it has "),
                    (prop_na*100), "% or more missing values: ",
                    paste0("\n - ", colnames(raw_data)[idx], collapse = ""))
      message(msg)
      raw_data <- raw_data[, !idx]
    }
  }
  
  # Check if there are columns with all NAs
  idx <- lapply(raw_data, function(x) all(is.na(x))) == TRUE
  # idx <- unname(idx[!(idx %in% excluded_columns)])
  if (sum(idx) > 0) {
    msg <- paste0("The following trait",
                  ifelse(sum(idx) > 1, "s were ", " was "),
                  "dropped because ",
                  ifelse(sum(idx) > 1, "they have ", "it has "),
                  "100% missing values: ",
                  paste0("\n - ", colnames(raw_data)[idx], collapse = ""))
    message(msg)
    raw_data <- raw_data[, !idx]
  }
  return(raw_data)
}

#' Assess normality of traits
#' 
#' @description 
#' Assess normality of traits in a data frame. 
#' 
#' @details 
#' The normality of each trait is assessed using a \emph{Shapiro-Wilk} test, 
#' under the following hypotheses:
#' 
#' \itemize{
#'   \item \eqn{H_0:} the sample comes from a normally distributed population.
#'   \item \eqn{H_1:} the sample does not come from a normally distributed 
#'   population.
#' }
#' 
#' Using a significance level of \eqn{\alpha = 0.05}. If the conclusion is that 
#' the sample does not come from a normally distributed population, then a 
#' number of transformations are performed, based on the transformation values 
#' passed with \code{transf_vals}. By default, the following transformation 
#' values are used \code{a = c(2, exp(1), 3, 4, 5, 6, 7, 8, 9, 10)} with the 
#' logarithmic (\code{log_a(x)}), power (\code{x^a}), and 
#' radical/root (\code{x^(1/a)}) functions.
#' 
#' @importFrom foreach %dopar%
#' @importFrom stats shapiro.test
#' 
#' @param raw_data Data frame containing the raw data.
#' @param excluded_columns Numeric vector containing the indices of the dataset 
#'     properties that are non-numeric, excluded columns.
#' @param cpus Number of CPUs to be used in the computation.
#' @param out_prefix Prefix for output files and plots.
#' @param plots_dir Path to the directory where plots should be stored.
#' @param transf_vals Numeric vector with the transformation values.
#' @param alpha Significance level.
#' @param pareto_scaling Boolean flag to indicate whether or not perform a 
#'     Pareto scaling on the normalised data.
#' @param show_stats Boolean flag to indicate whether or not to show the 
#'     normality assessment statistics (how many traits are normal, how many
#'     were transformed/normalised, and which transformations were applied).
#'     
#' @return List of data frames for the normal (\code{norm}) and skewed 
#' (\code{skew}) traits.
#' 
#' @export
#' 
#' @examples
#' \donttest{
#' # Toy dataset
#' example_data <- data.frame(ID = c(1,2,3,4,5), 
#'                            P1 = c("one", "two", "three", "four", "five"), 
#'                            T1 = rnorm(5), 
#'                            T2 = rnorm(5))
#' example_data_normalised <- MetaPipe::assess_normality(example_data, c(1, 2))
#' example_data_norm <- example_data_normalised$norm
#' example_data_skew <- example_data_normalised$skew
#' 
#' # Normal traits
#' knitr::kable(example_data_norm)
#' 
#' # Skewed traits (empty)
#' # knitr::kable(example_data_skew)
#' 
#' 
#' # F1 Seedling Ionomics dataset
#' data(ionomics) # Includes some missing data
#' ionomics_rev <- MetaPipe::replace_missing(ionomics, 
#'                                           excluded_columns = c(1, 2),
#'                                           replace_na =  TRUE)
#' ionomics_normalised <- 
#'   MetaPipe::assess_normality(ionomics_rev,
#'                              excluded_columns = c(1, 2),
#'                              out_prefix = "ionomics",
#'                              transf_vals = c(2, exp(1)))
#'                              
#' ionomics_norm <- ionomics_normalised$norm
#' ionomics_skew <- ionomics_normalised$skew
#' 
#' # Normal traits
#' knitr::kable(ionomics_norm[1:5, ])
#' 
#' # Skewed traits (partial output)
#' knitr::kable(ionomics_skew[1:5, 1:8])
#' 
#' # Clean up example outputs
#' MetaPipe:::tidy_up(c("HIST_", "ionomics_", "metapipe_"))
#' }
assess_normality <- function(raw_data, 
                             excluded_columns, 
                             cpus = 1, 
                             out_prefix = "metapipe", 
                             plots_dir = tempdir(), 
                             transf_vals = c(2, 
                                             exp(1), 
                                             3, 
                                             4, 
                                             5, 
                                             6, 
                                             7, 
                                             8, 
                                             9, 
                                             10),
                             alpha = 0.05,
                             pareto_scaling = FALSE,
                             show_stats = TRUE) {
  
  # Call core function: assess normality and transform traits
  raw_data_normalised <- assess_normality_core(raw_data, 
                                               excluded_columns,
                                               cpus, 
                                               out_prefix, 
                                               plots_dir, 
                                               transf_vals, 
                                               alpha)
  # Call the post-processing function
  raw_data_normalised_post <- 
    assess_normality_postprocessing(raw_data,
                                    excluded_columns,
                                    raw_data_normalised,
                                    out_prefix,
                                    pareto_scaling)
  
  # Show stats of the normalisation process
  if (show_stats) {
    assess_normality_stats(out_prefix)
  }
  
  return(raw_data_normalised_post)
}

#' QTL mapping
#' 
#' Perform QTL mapping using the \code{\link[qtl:scanone]{qtl:scanone(...)}} 
#' function to obtain LOD scores for all traits, peak positions, and markers.
#' 
#' @importFrom foreach %dopar%
#' 
#' @param x_data Cross-data frame containing genetic map data and traits.
#' @param cpus Number of CPUs to be used in the computation.
#' @param ... Arguments passed on to 
#'     \code{\link[qtl:scanone]{qtl::scanone}}.
# @inheritDotParams qtl::scanone -cross -pheno.col
#' 
#' @return Data frame containing the LOD scores.
#' @export
#'
#' @examples
#' \donttest{
#' # Create temp dir
#' tmp <- tempdir()
#' dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
#' 
#' # Toy dataset
#' excluded_columns <- c(1, 2)
#' population <- 5
#' seed <- 123
#' set.seed(seed)
#' example_data <- data.frame(ID = 1:population,
#'                            P1 = c("one", "two", "three", "four", "five"),
#'                            T1 = rnorm(population),
#'                            T2 = rnorm(population))
#' 
#' output <- MetaPipe::assess_normality(example_data, 
#'                                      excluded_columns, 
#'                                      show_stats = FALSE,
#'                                      out_prefix = paste0(tmp, "/tmp"))
#' 
#' # Create and store random genetic map (for testing only)
#' genetic_map <- MetaPipe:::random_map(population = population, 
#'                                      seed = seed)
#' # Load cross file with genetic map and raw data for normal traits
#' x <- MetaPipe::read.cross(genetic_map, output$norm)
#' 
#' x <- qtl::calc.genoprob(x, step = 1, error.prob = 0.001)
#' x_scone <- MetaPipe::qtl_scone(x, 1, model = "normal", method = "hk")
#' 
#' # F1 Seedling Ionomics dataset
#' data(ionomics) # Includes some missing data
#' data(father_riparia) # Genetic map
#' ionomics_rev <- MetaPipe::replace_missing(ionomics, 
#'                                           excluded_columns = c(1, 2),
#'                                           replace_na =  TRUE,
#'                                           out_prefix = paste0(tmp, "/tmp"))
#' ionomics_normalised <- 
#'   MetaPipe::assess_normality(ionomics_rev,
#'                              excluded_columns = c(1, 2),
#'                              out_prefix = file.path(tmp, "ionomics"),
#'                              transf_vals = c(2, exp(1)),
#'                              show_stats = FALSE)
#'                              
#' # Load cross file with genetic map and raw data for normal traits
#' x <- MetaPipe::read.cross(father_riparia, 
#'                           ionomics_normalised$norm,
#'                           genotypes = c("nn", "np", "--"))
#'                           
#' set.seed(seed)
#' x <- qtl::jittermap(x)
#' x <- qtl::calc.genoprob(x, step = 1, error.prob = 0.001)
#' 
#' x_scone <- MetaPipe::qtl_scone(x, 1, model = "normal", method = "hk")
#' 
#' # Clean temporal directory
#' # unlink(tmp, recursive = TRUE, force = TRUE)
#' MetaPipe:::tidy_up(tmp)
#' }
#' @family QTL mapping functions
qtl_scone <- function(x_data, cpus = 1, ...) {
  i <- NULL # Local binding
  # Start parallel backend
  cl <- parallel::makeCluster(cpus)# , setup_strategy = "sequential")
  doParallel::registerDoParallel(cl)
  
  # Load binary operator for backend
  # `%dopar%` <- foreach::`%dopar%`
  
  # Compute trait indices, accounting for the offset of ID and properties
  trait_indices <- 2:ncol(x_data$pheno)
  
  # Extract trait names
  traits <- colnames(x_data$pheno)
  
  x_scone <- foreach::foreach(i = trait_indices,
                     .combine = cbind) %dopar% {
                       # Run single scan
                       scone <- qtl::scanone(x_data, pheno.col = i, ...)
                       if(i == 2) {
                         record <- data.frame(
                           chr = scone$chr,
                           pos = scone$pos,
                           lod = scone$lod,
                           row.names = rownames(scone)
                         )
                         colnames(record)[3] <- traits[i]
                       }
                       else {
                         record <- data.frame(data = scone$lod)
                         colnames(record) <- traits[i]
                       }
                       record
                     }
  parallel::stopCluster(cl) # Stop cluster
  return(x_scone)
}

#' QTL mapping permutation test
#' 
#' Perform a QTL mapping permutation test using the 
#' \code{\link[qtl:scanone]{qtl:scanone(...)}} function to find significant QTL.
#' 
#' @importFrom foreach %dopar%
#' @importFrom graphics abline legend plot
#' @importFrom stats as.formula
#' 
#' @param x_data Cross-data frame containing genetic map data and traits.
#' @param cpus Number of CPUs to be used in the computation.
#' @param qtl_method QTL mapping method.
#' @param raw_data_normalised Normalised raw data, see 
#'     \code{\link{assess_normality}}.
#' @param lod_threshold LOD score threshold to look up for significant QTLs
#' @param parametric Boolean flag to indicate whether or not \code{x_data} 
#'     contains parametric (normal) traits.
#' @param n_perm Number of permutations.
#' @param plots_dir Output directory for plots.
#' @param ... Arguments passed on to 
#'     \code{\link[qtl:scanone]{qtl::scanone}}.
# @inheritDotParams qtl::scanone -cross -pheno.col -n.perm
#' 
#' @return Data frame containing the significant QTLs information.
#' @export
#' 
#' @examples
#' \donttest{
#' # Create temp dir
#' tmp <- tempdir()
#' dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
#' 
#' # Toy dataset
#' excluded_columns <- c(1, 2)
#' population <- 5
#' seed <- 123
#' set.seed(seed)
#' example_data <- data.frame(ID = 1:population,
#'                            P1 = c("one", "two", "three", "four", "five"),
#'                            T1 = rnorm(population),
#'                            T2 = rnorm(population))
#' 
#' output <- MetaPipe::assess_normality(example_data, 
#'                                      excluded_columns, 
#'                                      show_stats = FALSE,
#'                                      out_prefix = paste0(tmp, "/tmp"))
#' 
#' # Create and store random genetic map (for testing only)
#' genetic_map <- MetaPipe:::random_map(population = population, 
#'                                      seed = seed)
#' 
#' # Load cross file with genetic map and raw data for normal traits
#' x <- MetaPipe::read.cross(genetic_map, output$norm)
#' 
#' x <- qtl::calc.genoprob(x, step = 1, error.prob = 0.001)
#' x_scone <- MetaPipe::qtl_scone(x, 1, model = "normal", method = "hk")
#' x_qtl_perm <- MetaPipe::qtl_perm_test(x, 
#'                                       n_perm = 5, 
#'                                       model = "normal", 
#'                                       method = "hk",
#'                                       plots_dir = tmp)
#' x_qtl_perm_1000 <- MetaPipe::qtl_perm_test(x, 
#'                                            n_perm = 1000, 
#'                                            model = "normal", 
#'                                            method = "hk",
#'                                            plots_dir = tmp)
#' 
#' # F1 Seedling Ionomics dataset
#' data(ionomics) # Includes some missing data
#' data(father_riparia) # Genetic map
#' ionomics_rev <- MetaPipe::replace_missing(ionomics, 
#'                                           excluded_columns = c(1, 2),
#'                                           replace_na =  TRUE,
#'                                           out_prefix = paste0(tmp, "/tmp"))
#' ionomics_normalised <- 
#'   MetaPipe::assess_normality(ionomics_rev,
#'                              excluded_columns = c(1, 2),
#'                              out_prefix = file.path(tmp, "ionomics"),
#'                              transf_vals = c(2, exp(1)),
#'                              show_stats = FALSE)
#' 
#' # Load cross file with genetic map and raw data for normal traits
#' x <- MetaPipe::read.cross(father_riparia, 
#'                           ionomics_normalised$norm,
#'                           genotypes = c("nn", "np", "--"))
#'                           
#' set.seed(seed)
#' x <- qtl::jittermap(x)
#' x <- qtl::calc.genoprob(x, step = 1, error.prob = 0.001)
#' 
#' x_scone <- MetaPipe::qtl_scone(x, 1, model = "normal", method = "hk")
#' x_qtl_perm <- MetaPipe::qtl_perm_test(x, 
#'                                       n_perm = 5, 
#'                                       model = "normal", 
#'                                       method = "hk",
#'                                       plots_dir = tmp)
#' x_qtl_perm_1000 <- MetaPipe::qtl_perm_test(x, 
#'                                            n_perm = 1000, 
#'                                            model = "normal", 
#'                                            method = "hk",
#'                                            plots_dir = tmp)
#' 
#' # Clean temporal directory
#' # unlink(tmp, recursive = TRUE, force = TRUE)
#' MetaPipe:::tidy_up(tmp)
#' }
#' 
#' @family QTL mapping functions
qtl_perm_test <- function(x_data, 
                          cpus = 1, 
                          qtl_method = "par-scanone", 
                          raw_data_normalised = NULL, 
                          lod_threshold = 3, 
                          parametric = TRUE, 
                          n_perm = 1000, 
                          plots_dir = tempdir(), 
                          ...) {
  i <- NULL # Local bindings
  # Start parallel backend
  cl <- parallel::makeCluster(cpus)# , setup_strategy = "sequential")
  doParallel::registerDoParallel(cl)
  
  # Load binary operator for backend
  # `%dopar%` <- foreach::`%dopar%`
  
  # Compute trait indices, accounting for the offset of ID and properties
  trait_indices <- 2:ncol(x_data$pheno)
  
  # Extract trait names
  traits <- colnames(x_data$pheno)
  
  # Obtain number of individuals (population)
  num_indv <- nrow(x_data$pheno) #summary(x_data)[[2]]
  
  x_sum_map <- 
    foreach::foreach(i = trait_indices,
                     .combine = rbind) %dopar% {
                       if (!is.null(raw_data_normalised)) {
                         transf_info <- raw_data_normalised$trait == traits[i]
                         transf_info <- 
                           raw_data_normalised[transf_info, 
                                               c("transf", "transf_val")][1, ]
                       } else {
                         transf_info <- data.frame(transf = NA, transf_val = NA)
                       }
                       
                       # Structure for QTL
                       record <- data.frame(
                         # ID = i - 1,
                         qtl_ID = NA,
                         trait = traits[i],
                         ind = num_indv,
                         lg = NA,
                         lod_peak = NA,
                         pos_peak = NA,
                         marker = NA,
                         pos_p95_bay_int = NA,
                         marker_p95_bay_int = NA,
                         pvar = NA,
                         est_add = NA,
                         est_dom = NA,
                         p5_lod_thr = NA,
                         p10_lod_thr = NA,
                         pval = NA,
                         transf = transf_info$transf,
                         transf_val = transf_info$transf_val,
                         method = qtl_method,
                         p5_qtl = FALSE,
                         p10_qtl = FALSE
                       )
                       
                       # Run single scan
                       x_scone <- qtl::scanone(x_data, pheno.col = i, ...)
                       sum_x_scone <- summary(x_scone, 
                                              threshold = lod_threshold)
                       lod_cnt <- nrow(sum_x_scone)
                       if (!is.null(lod_cnt) && lod_cnt > 0) {
                         for (k in 1:lod_cnt) {
                           if (k > 1) {
                             # Create copy of record object
                             nrecord <- record[1, ]
                           } else {
                             # Copy record structured and data
                             nrecord <- record
                           }
                           
                           # Extract Peak QTL information
                           nrecord$lg <- sum_x_scone[k, "chr"]       
                           nrecord$lod_peak <- sum_x_scone[k, "lod"]
                           nrecord$pos_peak <- sum_x_scone[k, "pos"]
                           marker <- rownames(sum_x_scone)[k]
                           
                           # Verify if current QTL has a pseudomarker
                           marker_info <- 
                             transform_pseudo_marker(x_data, 
                                                     marker, 
                                                     nrecord$lg, 
                                                     nrecord$pos_peak)
                           nrecord$marker <- marker_info[1]
                           nrecord$pos_peak <- as.numeric(marker_info[2])
                           
                           # Create QTL ID: trait:LG@position
                           if (!is.na(nrecord$lg)) {
                             nrecord$qtl_ID <- 
                               with(nrecord, 
                                    sprintf("%s:%s@%f", 
                                            traits[i], 
                                            lg, 
                                            pos_peak))
                           }
                           
                           # Compute the 95% Bayes' CI
                           p95_bayes <- qtl::bayesint(x_scone, 
                                                      chr = nrecord$lg, 
                                                      expandtomarkers = TRUE, 
                                                      prob = 0.95)
                           p95_bayes <- unique(p95_bayes)
                           low_bound <- 1 #p95_bayes$pos == min(p95_bayes$pos)
                           upper_bound <- which.max(p95_bayes$pos) #p95_bayes$pos == max(p95_bayes$pos)
                           # Add new column for markers, 
                           # prevent duplicated row names
                           p95_bayes$marker <- NA 
                           
                           # Verify if the 95% Bayes' CI QTLs have pseudomarkers
                           for (l in 1:nrow(p95_bayes)) {
                             marker <- rownames(p95_bayes)[l]
                             marker_info <- 
                               transform_pseudo_marker(x_data, 
                                                       marker, 
                                                       p95_bayes[l, "chr"], 
                                                       p95_bayes[l, "pos"])
                             p95_bayes[l, "marker"] <- marker_info[1]
                             p95_bayes[l, "pos"] <- as.numeric(marker_info[2])
                           }
                           
                           nrecord$pos_p95_bay_int <- 
                             paste0(p95_bayes[low_bound, "pos"], "-",
                                    p95_bayes[upper_bound, "pos"])
                           nrecord$marker_p95_bay_int <- 
                             paste0(p95_bayes[low_bound, "marker"], "-",
                                    p95_bayes[upper_bound, "marker"])
                           
                           if (k > 1) {
                             record <- rbind(record, nrecord)
                           } else {
                             record <- nrecord
                           }
                         }
                         
                         x_scone_perm <- qtl::scanone(x_data, 
                                                      pheno.col = i, 
                                                      n.perm = n_perm, 
                                                      ...)
                         p5 <- summary(x_scone_perm)[[1]]  #  5% percent
                         p10 <- summary(x_scone_perm)[[2]] # 10% percent
                         
                         lod_plot <- save_plot(
                           {
                             plot(x_scone, ylab = "LOD Score")
                             abline(h = p5, 
                                    lwd = 2,
                                    lty = "solid",
                                    col = "red"
                             )
                             abline(h = p10, 
                                    lwd = 2, 
                                    lty = "dashed",
                                    col = "blue"
                             )
                             legend("topleft", 
                                    legend = c("5 %", "10 %"), 
                                    col = c("red", "blue"), 
                                    lty = c("solid", "dashed"), 
                                    lwd = 2
                                   )
                           },
                           paste0(plots_dir, "/LOD-", traits[i]),
                           width = 18
                         )
                         
                         record[, ]$p5_lod_thr <- p5
                         record[, ]$p10_lod_thr <- p10
                         
                         p5_idx <- which(record$lod_peak >= p5)
                         p10_idx <- which(record$lod_peak >= p10)
                         if (length(p5_idx) > 0) {
                           record[p5_idx, ]$p5_qtl <- TRUE
                         }
                         if (length(p10_idx) > 0) {
                           record[p10_idx, ]$p10_qtl <- TRUE
                         }
                         
                         # For parametric QTL mapping only
                         if (parametric) {
                           chr <- as.numeric(sum_x_scone$chr)
                           pos <- as.numeric(sum_x_scone$pos)
                           qtl_s <- qtl::makeqtl(x_data, 
                                                 chr, 
                                                 pos, 
                                                 what = c("prob"))
                           
                           for (m in 1:length(chr)) {
                             #qtl_s <- makeqtl(x_data, chr[m], pos[m], what=c("prob"))
                             #f <- as.formula(paste0("y~",paste0("Q",seq(1:nrow(sum_x_scone)), collapse = " + ")))
                             f <- as.formula(paste0("y~", 
                                                    paste0("Q", 
                                                           m, 
                                                           collapse = " + ")))
                             fit_qtl <- qtl::fitqtl(x_data, 
                                                    pheno.col = i, 
                                                    qtl_s, 
                                                    formula = f, 
                                                    get.ests = TRUE, 
                                                    ...)
                             sum_fit_qtl <- summary(fit_qtl)
                             
                             if (length(sum_fit_qtl) > 0) {
                               p.var <- as.numeric(sum_fit_qtl[[1]][1, "%var"])
                               pvalue.f <- 
                                 as.numeric(sum_fit_qtl[[1]][, "Pvalue(F)"])[1]
                               estimates <- 
                                 as.numeric(sum_fit_qtl$ests[, "est"])[-1]
                               record[m, ]$pvar <- p.var
                               record[m, ]$pval <- pvalue.f
                               record[m, ]$est_add <- estimates[1]
                               record[m, ]$est_dom <- estimates[2]
                             }
                           }
                         }
                       }
                       record
                     }
  parallel::stopCluster(cl) # Stop cluster
  return(x_sum_map)
}
