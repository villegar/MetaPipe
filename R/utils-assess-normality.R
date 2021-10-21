#' Assess normality of traits
#' 
#' @description 
#' Assess normality of traits in a data frame. 
#' 
#' @details 
#' The normality of each trait is checked using a \emph{Shapiro-Wilk} test, 
#' under the following hypotheses:
#' 
#' \itemize{
#'   \item \eqn{H_0:} the sample comes from a normally distributed population
#'   \item \eqn{H_1:} the sample does not come from a normally distributed 
#'   population
#' }
#' 
#' Using a significance level of \eqn{\alpha = 0.05}. If the conclusion is that 
#' the sample does not come from a normally distributed population, then a 
#' number of transformations are performed, based on the transformation values 
#' passed with \code{transf_vals}. By default, the following transformation 
#' values are used \code{c(2, exp(1), 3, 4, 5, 6, 7, 8, 9, 10)} with the 
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
#'
#' @return Structure containing the normalised data, if a suitable 
#' transformation was found, otherwise returns the original data.
#'
#' @examples
#' \donttest{
#' # Toy dataset
#' example_data <- data.frame(ID = c(1,2,3,4,5), 
#'                            P1 = c("one", "two", "three", "four", "five"), 
#'                            T1 = rnorm(5), 
#'                            T2 = rnorm(5))
#' example_data_normalised <- MetaPipe:::assess_normality_core(example_data, 
#'                                                             c(1, 2))
#' knitr::kable(example_data_normalised)
#' 
#' 
#' # F1 Seedling Ionomics dataset
#' data(ionomics) # Includes some missing data
#' ionomics_rev <- MetaPipe::replace_missing(ionomics, 
#'                                           excluded_columns = c(1, 2),
#'                                           replace_na =  TRUE)
#' ionomics_normalised <- 
#'   MetaPipe:::assess_normality_core(ionomics_rev,
#'                                    excluded_columns = c(1, 2),
#'                                    transf_vals = c(2, exp(1)))
#' # Show one entry for each of the first ten traits (left to right)
#' knitr::kable(ionomics_normalised[nrow(ionomics) * c(1:10), ])
#' }
#' 
#' @seealso \code{\link{assess_normality_postprocessing}} and 
#' \code{\link{assess_normality_stats}}
#' 
#' @keywords internal
#' @noRd
assess_normality_core <- function(raw_data, 
                             excluded_columns, 
                             cpus = 1, 
                             out_prefix = file.path(tempdir(), "metapipe"), 
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
                             alpha = 0.05) {
  i <- NULL # Local binding
  # Start parallel backend
  cl <- parallel::makeCluster(cpus)# , setup_strategy = "sequential")
  doParallel::registerDoParallel(cl)
  
  # Load binary operator for backend
  #`%dopar%` <- foreach::`%dopar%`
  
  # Verify that plots_dir exists
  if (!dir.exists(plots_dir))
    dir.create(plots_dir, FALSE)
  
  # Exclude column 1, ID
  excluded_columns <- unique(c(1, excluded_columns))
  
  # Ignore ID and properties
  raw_data <- raw_data[, -excluded_columns]
  
  # Compute trait indices, accounting for the offset of ID and properties
  trait_indices <- 1:ncol(raw_data) 
  
  # Extract traits (column names)
  traits <- colnames(raw_data)
  raw_data_normalised <- 
    foreach::foreach(i = trait_indices, 
                     .combine = rbind) %dopar% {
                       # Create en empty entry for current trait
                       record <- data.frame(
                         index = i,
                         trait = traits[i],
                         values = raw_data[, i],
                         flag = "Skewed",
                         transf = "",
                         transf_val = NA,
                         stringsAsFactors = FALSE
                       )
                       
                       # Verify the current trait has at least 3 non-NA rows
                       if (sum(is.finite(raw_data[, i]), na.rm = TRUE) > 2) {
                         # Assess normality of trait before transforming it
                         pvalue <- shapiro.test(raw_data[, i])[[2]]
                         if(pvalue <= alpha) { # Data must be transformed
                           tmp <- MetaPipe::transform_data(
                             data = raw_data[, i],
                             trait = traits[i],
                             alpha = alpha,
                             index = i,
                             transf_vals = transf_vals,
                             plots_prefix = paste0(plots_dir, "/HIST")
                           )
                           
                           if (length(tmp)) {
                             tmp$index <- i
                             tmp$flag <- "Normal"
                             record <- tmp
                           }
                         } else { # Normal data
                           xlab <- traits[i]
                           transformation <- "NORM"
                           prefix <- paste0(plots_dir, 
                                            "/HIST_", 
                                            i, 
                                            "_", 
                                            transformation)
                           generate_hist(data = raw_data[, i], 
                                         title = traits[i], 
                                         prefix = prefix, 
                                         xlab = xlab, 
                                         is_trait = TRUE)
                           record$flag <- "Normal"
                         }
                       }
                       record
                     }
  
  parallel::stopCluster(cl) # Stop cluster
  return(raw_data_normalised)
}

#' Post-processing for \code{\link{assess_normality_core}}
#' 
#' Post-processing for the normality assessment of the traits.
#' 
#' It creates four CSV files containing the following data:
#' \itemize{
#'   \item \code{paste0(out_prefix, "_raw_data_normalised_all.csv")}: 
#'   data frame in long format created with the \code{\link{assess_normality_core}}
#'   function. Contains all the traits, including name, values, and 
#'   transformation (if applicable).
#'   \item \code{paste0(out_prefix, "_raw_data_norm.csv")}: data frame in wide
#'   format containing traits with a normal distribution.
#'   \item \code{paste0(out_prefix, "_raw_data_non_par.csv")}: data frame in 
#'   wide format containing traits \emph{without} a normal distribution.
#'   \item \code{paste0(out_prefix, "_normalisation_stats.csv")}: data frame
#'   containing details of the \code{\link{assess_normality_core}} process, like 
#'   the total number of traits transformed and normalisations performed.
#' }
#' 
#' @param raw_data Data frame containing the raw data.
#' @param excluded_columns Numeric vector containing the indices of the dataset 
#'     properties that are non-numeric, excluded columns.
#' @param raw_data_normalised Data frame containing the normalised raw data, 
#'     created with \code{\link{assess_normality_core}}.
#' @param out_prefix Prefix for output files and plots.
#' @param pareto_scaling Boolean flag to indicate whether or not perform a 
#'     Pareto scaling on the normalised data.
#'
#' @return List of data frames for the normal (\code{norm}) and skewed 
#' (\code{skew}) traits.
#'
#' @examples
#' \donttest{
#' # Toy dataset
#' example_data <- data.frame(ID = c(1,2,3,4,5), 
#'                            P1 = c("one", "two", "three", "four", "five"), 
#'                            T1 = rnorm(5), 
#'                            T2 = rnorm(5))
#' example_data_normalised <- MetaPipe:::assess_normality_core(example_data, 
#'                                                             c(1, 2))
#' example_data_normalised_post <- 
#'  MetaPipe:::assess_normality_postprocessing(example_data, 
#'                                             c(1, 2), 
#'                                             example_data_normalised)
#' example_data_norm <- example_data_normalised_post$norm
#' example_data_skew <- example_data_normalised_post$skew
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
#'   MetaPipe:::assess_normality_core(ionomics_rev,
#'                                    excluded_columns = c(1, 2),
#'                                    transf_vals = c(2, exp(1)))
#' ionomics_normalised_post <- 
#'   MetaPipe:::assess_normality_postprocessing(ionomics_rev, 
#'                                              c(1, 2), 
#'                                              ionomics_normalised)
#' ionomics_norm <- ionomics_normalised_post$norm
#' ionomics_skew <- ionomics_normalised_post$skew
#' # Normal traits
#' knitr::kable(ionomics_norm[1:5, ])
#' 
#' # Skewed traits (partial output)
#' knitr::kable(ionomics_skew[1:5, 1:10])
#' }
#' 
#' @seealso \code{\link{assess_normality_core}} and 
#' \code{\link{assess_normality_stats}}
#' 
#' @keywords internal
#' @noRd
assess_normality_postprocessing <- function(raw_data, 
                                            excluded_columns,
                                            raw_data_normalised,
                                            out_prefix = file.path(tempdir(), 
                                                                   "metapipe"), 
                                            pareto_scaling = FALSE) {
  trait <- transf_val <- NULL # Local binding
  # Verify the raw_data_normalised object has the right structure
  expected_columns <- c("index", 
                        "trait", 
                        "values", 
                        "flag", 
                        "transf", 
                        "transf_val")
  if (!all(expected_columns %in% colnames(raw_data_normalised)))
    stop("raw_data must be the output of the function assess_normality_core")
  
  # Exclude column 1, ID
  excluded_columns <- unique(c(1, excluded_columns))
  
  # Separate normal and non-parametric entries from the normalised data
  raw_data_normalised_norm <- 
    raw_data_normalised[raw_data_normalised$flag == "Normal", ]
  raw_data_normalised_non_par <- 
    raw_data_normalised[raw_data_normalised$flag == "Skewed", ]
  
  # Extract trait names for both normal and non-parametric data
  traits_non_par <- unique(as.character(raw_data_normalised_non_par$trait))
  traits_norm <- unique(as.character(raw_data_normalised_norm$trait))
  traits_non_par_len <- length(traits_non_par)
  traits_norm_len <- length(traits_norm)
  
  # Create new objects with the normalised data for normal traits and original
  # raw data for non-parametric traits
  raw_data_non_par <- NULL
  raw_data_norm <- NULL
  if (traits_non_par_len > 0)
    raw_data_non_par <- raw_data[, traits_non_par]
  if (traits_norm_len > 0 ) {
    raw_data_norm <- 
      data.frame(matrix(vector(), 
                        nrow(raw_data_normalised_norm) / traits_norm_len, 
                        traits_norm_len,
                        dimnames = list(c(), traits_norm)),
                 stringsAsFactors = FALSE)
    
    for (i in 1:traits_norm_len) {
      raw_data_norm[i] <- 
        subset(raw_data_normalised_norm, trait == traits_norm[i])$values
    }
  }
  
  # Append excluded columns for scaling, if pareto_scaling == TRUE
  ## Normal traits
  if (!is.null(raw_data_norm)) {
    raw_data_norm <-
      cbind(raw_data[, 1, drop = FALSE],
            if (pareto_scaling)
              paretoscale(raw_data_norm)
            else
              raw_data_norm)
  }
  
  ## Skewed traits
  if (!is.null(raw_data_non_par)) {
    raw_data_non_par <-
      cbind(raw_data[, 1, drop = FALSE],
            if (pareto_scaling)
              paretoscale(raw_data_non_par)
            else
              raw_data_non_par
      )
  }
  
  # Generate basic stats from the normalisation process
  raw_data_rows <- nrow(raw_data)
  norm_traits_normalised_count <-
    nrow(raw_data_normalised_norm[raw_data_normalised_norm$transf != "",])
  norm_traits_normalised_count <- norm_traits_normalised_count / raw_data_rows
  norm_traits_count <- nrow(raw_data_normalised_norm) / raw_data_rows
  total_traits <- nrow(raw_data_normalised)/raw_data_rows
  transformations <- unique(raw_data_normalised[c("transf", "transf_val")])
  # Drop blank transformation, NULL transformation
  # transformations <- transformations[-nrow(transformations), ]
  sorting <- order(transformations$transf, transformations$transf_val)
  transformations <- transformations[sorting, ]
  
  # Create data frame containing the stats
  normalisation_stats <- 
    data.frame(key = c("total", "norm_traits", "norm_traits_normalised"),
               values = c(total_traits, 
                          norm_traits_count, 
                          norm_traits_normalised_count),
               stringsAsFactors = FALSE)
  if (nrow(transformations) > 0) {
    for (i in 1:nrow(transformations)) {
      if(is.na(transformations$transf_val[i]))
        next
      key <-
        paste0(transformations$transf[i], "\t", transformations$transf_val[i])
      tmp <- 
        subset(raw_data_normalised_norm,
               raw_data_normalised_norm$transf == transformations$transf[i])
      tmp <- subset(tmp, transf_val == transformations$transf_val[i])
      value <- nrow(tmp) / raw_data_rows
      normalisation_stats <- rbind(normalisation_stats, c(key, value))
    }
  }
  
  # Write to disk new data structures
  write.csv(raw_data_normalised,
            file = paste0(out_prefix, "_raw_data_normalised_all.csv"),
            row.names = FALSE)
  write.csv(raw_data_norm,
            file = paste0(out_prefix, "_raw_data_norm.csv"),
            row.names = FALSE)
  write.csv(raw_data_non_par,
            file = paste0(out_prefix, "_raw_data_non_par.csv"),
            row.names = FALSE)
  write.csv(normalisation_stats,
            file = paste0(out_prefix, "_normalisation_stats.csv"),
            row.names = FALSE)
  
  return(list(norm = raw_data_norm, skew = raw_data_non_par))
}

#' Assess normality statistics
#' 
#' Statistics for the normality assessment of the traits. Uses an output file 
#' generated by \code{\link{assess_normality_postprocessing}}, therefore, it 
#' cannot be executed without previously running the post-processing function.
#'
#' @param out_prefix Prefix for output files and plots.
#'
#' @examples
#' \donttest{
#' # Toy dataset
#' example_data <- data.frame(ID = c(1,2,3,4,5), 
#'                            P1 = c("one", "two", "three", "four", "five"), 
#'                            T1 = rnorm(5), 
#'                            T2 = rnorm(5))
#' example_data_normalised <- MetaPipe:::assess_normality_core(example_data, 
#'                                                             c(1, 2))
#' example_data_normalised_post <- 
#'  MetaPipe:::assess_normality_postprocessing(example_data, 
#'                                            c(1, 2), 
#'                                            example_data_normalised)
#' MetaPipe:::assess_normality_stats()
#' 
#' 
#' # F1 Seedling Ionomics dataset
#' data(ionomics) # Includes some missing data
#' ionomics_rev <- MetaPipe::replace_missing(ionomics, 
#'                                           excluded_columns = c(1, 2),
#'                                           replace_na =  TRUE)
#' ionomics_normalised <- 
#'   MetaPipe:::assess_normality_core(ionomics_rev,
#'                                    excluded_columns = c(1, 2),
#'                                    transf_vals = c(2, exp(1)))
#' ionomics_normalised_post <- 
#'   MetaPipe:::assess_normality_postprocessing(ionomics_rev, 
#'                                             c(1, 2), 
#'                                             ionomics_normalised)
#' MetaPipe:::assess_normality_stats()
#' }
#' 
#' @seealso \code{\link{assess_normality_core}} and 
#' \code{\link{assess_normality_postprocessing}}
#' 
#' @keywords internal
#' @noRd
assess_normality_stats <- function(out_prefix = file.path(tempdir(), 
                                                          "metapipe")) {
  stats_filename <- paste0(out_prefix, "_normalisation_stats.csv")
  if (!file.exists(stats_filename))
    stop(paste0("The file ", stats_filename, " was not found"))
  
  # Loading stats for the normality assessment process
  normalisation_stats <- read.csv(stats_filename)
  total_traits <- normalisation_stats[1, 2]
  norm_traits_count <- normalisation_stats[2, 2]
  norm_traits_normalised_count <- normalisation_stats[3, 2]
  transformations <- normalisation_stats[-c(1:3),]
  
  # Create message for user with the summary
  msg <- paste0(sprintf("%-45s", "Total traits (excluding all NAs traits):"), 
                total_traits, "\n",
                sprintf("%-45s", "Normal traits (without transformation):"),
                (norm_traits_count - norm_traits_normalised_count), "\n",
                sprintf("%-45s", "Normal traits (transformed):"),
                norm_traits_normalised_count, "\n",
                sprintf("%-45s", "Total normal traits:"), 
                norm_traits_count, "\n",
                sprintf("%-45s", "Total skewed traits: "),
                (total_traits - norm_traits_count))
  
  if (nrow(transformations) > 0) {
    msg <- paste0(msg, "\n\nTransformations summary:")
    msg <- paste0(msg, 
                  sprintf("\n\t%-10s%-10s%-10s", "f(x)", "Value", "# traits"))
                  #"\n\tf(x)\tValue \t# traits")
    for (i in 1:nrow(transformations)) {
      tmp <- strsplit(as.character(transformations[i, 1]), '\t')[[1]]
      msg <- paste0(msg,
                    sprintf("\n\t%-10s%-10s%-10s", 
                            tmp[1], 
                            tmp[2], 
                            transformations[i, 2]))
                    # tmp[1], "\t", tmp[2],"\t", transformations[i, 2])
    }
  }
  msg <- paste0(msg, "\n") # Clean output
  message(msg)
}