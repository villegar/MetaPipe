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
