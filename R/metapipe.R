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
#' load_raw("example_data.csv", c(1, 2))
#' load_raw("example_data_dup.csv", c(1, 2))
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
#' example_data <- data.frame(ID = c(1,2,3,4,5), 
#'                            P1 = c("one", "two", "three", "four", "five"), 
#'                            F1 = rnorm(5), 
#'                            F2 = rnorm(5))
#' example_data$F1[2:3] <- NA 
#' example_data$F2[4] <- NA 
#' replace_missing(example_data, c(1, 2))
#' replace_missing(example_data, c(1, 2), prop_na =  0.25)
#' replace_missing(example_data, c(1, 2), replace_na =  TRUE)
replace_missing <- function(raw_data, excluded_columns, out_prefix = "metapipe", prop_na = 0.5, replace_na = FALSE) {
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
      msg <- "The following features were dropped because they have "
      msg <- paste0(msg, (prop_na*100), "% or more missing values: ")
      msg <- paste0(msg, paste(colnames(raw_data)[NACount], collapse = ", "), "\n")
      message(msg)
      raw_data[, NACount] <- NULL
    }
  }
  return(raw_data)
}

#' Assess normality of features
#' @importFrom foreach %dopar%
#' @param raw_data data frame containing the raw data
#' @param excluded_columns vector containing the indices of the data set properties, excluded columns
#' @param cpus number of CPUS to be used
#' @param out_prefix prefix for output files, plots, etc.
#' @param plots_dir path to the directory where plots should be stored
#' @param transf_vals transformation values
#'
#' @return Structure containing the normalised data, if a suitable transformation was found
#' @export
#'
#' @examples
#' example_data <- data.frame(ID = c(1,2,3,4,5), 
#'                            P1 = c("one", "two", "three", "four", "five"), 
#'                            F1 = rnorm(5), 
#'                            F2 = rnorm(5))
#' assess_normality(example_data, c(1, 2))
assess_normality <- function(raw_data, 
                             excluded_columns, 
                             cpus = 1, 
                             out_prefix = "metapipe", 
                             plots_dir = getwd(), 
                             transf_vals = c(2, exp(1), 3, 4, 5, 6, 7, 8, 9, 10)) {
  # Start parallel backend
  cl <- parallel::makeCluster(cpus)
  doParallel::registerDoParallel(cl)
  
  # Load binary operator for backend
  #`%dopar%` <- foreach::`%dopar%`
  
  # Exclude column 1, ID
  excluded_columns <- unique(c(1, excluded_columns))
  
  # Ignore ID and properties
  raw_data <- raw_data[, -excluded_columns]
  #len_excluded_columns <- length(excluded_columns)
  # Compute feature indices, accounting for the offset of ID and properties
  #feature_indices <- 1:(ncol(raw_data) - len_excluded_columns)
  feature_indices <- 1:ncol(raw_data) 
  # Extract features (column names)
  #features <- colnames(raw_data)[len_excluded_columns + feature_indices]
  features <- colnames(raw_data)
  raw_data_transformed <- foreach::foreach(i = feature_indices,
                                  .combine = rbind) %dopar% {
                                    # Create and populate entry for current feature
                                    record <- data.frame( 
                                      index = i,
                                      feature = features[i],
                                      values = raw_data[, i],
                                      flag = "Non-normal",
                                      transf = "",
                                      transf.value = NA
                                    )
                                    
                                    # Verify the current feature has at least 3 non-NA rows
                                    if(sum(is.finite(raw_data[, i]), na.rm = TRUE) > 2) {
                                      # Assess normality of feature before transforming it
                                      pvalue <- shapiro.test(raw_data[, i])[[2]]
                                      if(pvalue <= 0.05) { # Data must be transformed
                                        tmp <- transform_data(data = raw_data[, i], 
                                                              feature = features[i], 
                                                              alpha = 0.05, 
                                                              index = i, 
                                                              transf_vals = transf_vals, 
                                                              plots_prefix = paste0(plots_dir, "/HIST")
                                        )
                                        
                                        if(length(tmp)) {
                                          tmp$flag <- "Normal"
                                          record <- tmp
                                        }
                                      }
                                      else { # Normal data
                                        xlab <- features[i]
                                        transformation <- "NORM"
                                        prefix <- paste0(plots_dir,"/HIST_", i, "_", transformation)
                                        MetaPipe::generate_hist(data = raw_data[, i], 
                                                                feature = features[i], 
                                                                prefix = prefix, 
                                                                xlab = xlab)
                                        record$flag <- "Normal"
                                      }
                                    }
                                    record
                                  }
  
  parallel::stopCluster(cl) # Stop cluster
  raw_data_transformed$flag <- as.factor(raw_data_transformed$flag)
  return(raw_data_transformed)
}

