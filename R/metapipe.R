#' Load raw data
#' 
#' @param raw_data_filename filename containing the raw data, with or without full path
#' @param excluded_columns vector containing the indices of the data set properties, excluded columns
#'
#' @return pre-processed raw data [data frame]
#' @export
#'
#' @examples
#' example_data <- data.frame(ID = c(1,2,3,4,5), 
#'                            P1 = c("one", "two", "three", "four", "five"), 
#'                            F1 = rnorm(5), 
#'                            F2 = rnorm(5))
#' write.csv(example_data, "example_data.csv", row.names = FALSE)
#' write.csv(example_data[c(1:5, 1, 2), ], "example_data_dup.csv", row.names = FALSE)
#' load_raw("example_data.csv", c(2))
#' load_raw("example_data_dup.csv", c(2))
load_raw <- function(raw_data_filename, excluded_columns) {
  # Load and clean raw data
  raw_data <- read.csv(raw_data_filename)
  # Exclude column 1, ID
  excluded_columns <- unique(c(1, excluded_columns))
  
  # Aggregate data by row, computing the mean
  mean_raw_data <- aggregate(raw_data[, -excluded_columns],
                             by = list(raw_data[, 1]),
                             mean, 
                             na.action = na.omit)
  colnames(mean_raw_data)[1] <- colnames(raw_data)[1]
  
  # Joins the agreggated data with the excluded columns (properties)
  mean_raw_data <- dplyr::left_join(raw_data[, excluded_columns],
                                    mean_raw_data,
                                    by = colnames(raw_data)[1])
  
  # Find duplicated row values on the first column (ID)
  mean_raw_data <- mean_raw_data[!duplicated(mean_raw_data[, 1]), ]
  rownames(mean_raw_data) <- 1:nrow(mean_raw_data)
  return(mean_raw_data)
}

#' Replace missing values (NAs)
#'
#' @param raw_data data frame containing the raw data
#' @param excluded_columns vector containing the indices of the data set properties, excluded columns
#' @param out_prefix prefix for output files, plots, etc.
#' @param prop_na proportion of missing values, if a feature exceeds this threshold, then it is dropped out
#' @param replace_na boolean flag to indicate whether or not missing values should be replaced by half of the minimum value
#'
#' @return data frame containing the raw data without missing values
#' @export
#'
#' @examples
#' #' example_data <- data.frame(ID = c(1,2,3,4,5), 
#'                            P1 = c("one", "two", "three", "four", "five"), 
#'                            F1 = rnorm(5), 
#'                            F2 = rnorm(5))
#' example_data$F1[2:3] <- NA 
#' example_data$F2[4] <- NA 
#' replace_missing(example_data, c(2))
#' replace_missing(example_data, c(2), prop_na =  0.25)
#' replace_missing(example_data, c(2), replace_na =  TRUE)
replace_missing <- function(raw_data, excluded_columns, out_prefix = getwd(), prop_na = 0.5, replace_na = FALSE) {
  # Exclude column 1, ID
  excluded_columns <- unique(c(1, excluded_columns))
  
  # Missing values are replaced by half of the minimum non-zero value for each feature.
  if(replace_na) {
    NA2halfmin <- function(x) suppressWarnings(replace(x, is.na(x), (min(x, na.rm = TRUE)/2)))
    raw_data[,-excluded_columns] <- lapply(raw_data[,-excluded_columns], NA2halfmin)
  } else {
    NACount <- which(colMeans(is.na(raw_data[,-excluded_columns])) >= prop_na) + length(excluded_columns)
    if(length(NACount)) {
      write.csv(raw_data[, c(excluded_columns, NACount)], file = paste0(out_prefix,"_NA_raw_data.csv"), row.names = FALSE)
      cat(paste0("The following features were dropped because they have ", (prop_na*100), "% or more missing values: "))
      cat(paste0(paste(colnames(raw_data)[NACount], collapse = ", "), "\n"))
      raw_data[, NACount] <- NULL
    }
  }
  return(raw_data)
}
