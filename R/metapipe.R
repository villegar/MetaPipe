
load_raw <- function(raw_data_filename, excluded_columns) {
  # Load and clean raw data
  raw_data <- read.csv(raw_data) 
  ncols <- ncol(raw_data)
  len_excluded_columns <- length(excluded_columns)
  mean_raw_data <- aggregate(raw_data[, (len_excluded_columns + 1):ncols],
                             by = list(raw_data[, 1]),
                             mean, 
                             na.action = na.omit)
  colnames(mean_raw_data)[1] <- colnames(raw_data)[1]
  mean_raw_data <- left_join(raw_data[, excluded_columns],
                             mean_raw_data,
                             by = colnames(raw_data)[1])
  mean_raw_data <- mean_raw_data[!duplicated(mean_raw_data[, 1]),]
  rownames(mean_raw_data) <- 1:nrow(mean_raw_data)
  mean_raw_data_rows <- nrow(mean_raw_data)
}