#' Check variable types
#'
#' Check variable types in a data frame. Useful whenever the data must be 
#' filtered prior some numerical analysis
#' 
#' @param raw_data data frame containing the raw data 
#' @param excluded_columns list of excluded columns (won't be check)
#' @param numeric boolean to indicate the operation mode. If \code{TRUE} will
#'     look up for non-numeric variables and exclude them, otherwise it will
#'     filter numeric variables, including those with the following classes:
#'     \code{complex} (e.g. \code{1+i})
#'     \code{integer} (e.g. \code{1L})
#'     \code{numeric} (e.g. \code{1}, \code{1.01})
#' @param quiet boolean to hide warning messages
#'     
#' @return vector with the column indices that don't meet the checking criteria, 
#' only numeric (\code{numeric = TRUE}) or non-numeric (\code{numeric = FALSE}),  
#' if \code{excluded_columns != NULL}, returns a vector with both originally 
#' excluded and newly found column indices
#'
#' @examples
#' example_data <- data.frame(ID = c(1,2,3,4,5), 
#'                            P1 = c("one", "two", "three", "four", "five"), 
#'                            T1 = rnorm(5))
#' MetaPipe:::check_types(example_data)
#' MetaPipe:::check_types(example_data, numeric = FALSE)
#' MetaPipe:::check_types(example_data, excluded_columns = 1, numeric = FALSE)
#' MetaPipe:::check_types(example_data, quiet = FALSE)
#' 
#' @keywords internal
check_types <- function(raw_data,
                        excluded_columns = NULL,
                        numeric = TRUE,
                        quiet = TRUE) {
  # Extract data type (class) of each variable (column)
  var_types <- sapply(raw_data, typeof)
  if (!is.null(excluded_columns)) {
    var_types <- var_types[-excluded_columns]
  }
  if (numeric) { # if numeric = TRUE, find non-numeric columns
    idx <- !(var_types %in% c("numeric", "double", "integer", "complex"))
  } else { # if numeric = FALSE, find numeric columns
    idx <- var_types %in% c("numeric", "double", "integer", "complex")
  }
  # If any variables were found to meet the search criteria
  if (any(idx)) {
    if (!quiet) {
      warning(paste0("\nThe following variable",
                     ifelse(sum(idx) > 1, "s ", " "), # Plural
                     "will be excluded, as ",
                     ifelse(sum(idx) > 1, "they are ", "it is "), # Plural
                     ifelse(numeric, "non-", ""), # Check mode
                     "numeric: ",
                     paste0("\n - ", names(var_types)[idx], collapse = "")))
    }
    # Ammend indices to entire dataset
    idx <- colnames(raw_data) %in% names(var_types)[idx]
    # Append to the original list
    excluded_columns <- unique(c(excluded_columns, which(idx)))
  }
  return(excluded_columns)
}

#' Generate a pseudo-random genotypes
#' 
#' Generate a pseudo-random list of genotypes. Particularly useful for testing.
#'
#' @param genotypes Character vector containing the base genotypes.
#' @param size Output length.
#' @param seed Seed for reproducibility.
#'
#' @return Character vector containing the pseudo-random genotypes.
#'
#' @examples
#' MetaPipe:::random_genotypes()
#' MetaPipe:::random_genotypes(genotypes = c("nn", "np"))
#' MetaPipe:::random_genotypes(size = 3)
#' 
#' @seealso \code{\link{random_map}}
#' 
#' @keywords internal
random_genotypes <- function(genotypes = c("A", "H", "B"),
                             size = 100,
                             seed = NULL) {
  if (!is.null(seed))
    set.seed(seed)
  return(genotypes[sample(1:length(genotypes), size, replace = TRUE)])
}

#' Generate a pseudo-random genetic map
#' 
#' Generate a pseudo-random genetic map. Particularly useful for testing.
#'
#' @param genotypes Character vector containing the base genotypes.
#' @param lg Numeric vector containing the linkage groups.
#' @param markers Number of markers per linkage group.
#' @param population Population size (rows).
#' @param seed Seed for reproducibility.
#'
#' @return Data frame containing the pseudo-random genetic map.
#'
#' @examples
#' gmap_1 <- MetaPipe:::random_map()
#' gmap_2 <- MetaPipe:::random_map(genotypes = c("nn", "np"))
#' gmap_3 <- MetaPipe:::random_map(population = 3)
#' 
#' @seealso \code{\link{random_genotypes}}
#'  
#' @keywords internal
random_map <- function(genotypes = c("A", "H", "B"),
                       lg = 1:10,
                       markers = 10,
                       population = 100,
                       seed = NULL) {
  marker_names <- paste0(rep(paste0("S", lg, "_"), each = markers), 1:markers)
  # Temporal vector for LG and positions
  tmp <- data.frame(lg = rep(lg, each = markers), 
                    pos = rep(1:markers))
  
  # Empty map
  map <- data.frame(ID = c("", "", 1:population))
  for (k in 1:length(marker_names)) {
    if (!is.null(seed))
      seed <- seed + k
    new_genotypes <- c(tmp[k, 1],
                       tmp[k, 2],
                       random_genotypes(genotypes, population, seed))
    map <- cbind(map, new_genotypes)
  }
  
  # Add marker names (columns)
  colnames(map) <- c("ID", marker_names)
  
  return(map)
}

#' Replace NAs
#' 
#' Replace missing values (\code{NA}s) by half of the minimum value. If all the 
#' data points are missing, then returns the original input
#'
#' @param x original data points containing missing values (\code{NA}s)
#'
#' @return new data without missing entries, unless all of the original data
#' points were \code{NA}
#'
#' @examples
#' MetaPipe:::rplc_na(NA)
#' MetaPipe:::rplc_na(c(1, NA, 3, NA))
#' 
#' @keywords internal
rplc_na <- function(x) {
  replace(x, 
          is.na(x), 
          ifelse(all(is.na(x)), NA, (min(x, na.rm = TRUE) / 2))
  )
}

#' List files matching a pattern and delete them
#'
#' @param patterns Vector of strings with patterns to be deleted.
#' @param path Path to the directory containing the files
#' 
#' @keywords internal
tidy_up <- function(patterns, path = here::here()) {
  # List files
  filenames <- lapply(patterns, 
                      function(x) list.files(path, x, full.names = TRUE))
  suppressWarnings(
    . <- lapply(filenames, file.remove)
  )
}

