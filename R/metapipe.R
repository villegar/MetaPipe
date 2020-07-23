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
  raw_data <- read.csv(raw_data_filename, stringsAsFactors = FALSE)
  # Exclude column 1, ID
  excluded_columns <- unique(c(1, excluded_columns))
  
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
#' 
#' @seealso \code{\link{assess_normality_postprocessing}} and \code{\link{assess_normality_stats}}
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
  `%dopar%` <- foreach::`%dopar%`
  
  # Exclude column 1, ID
  excluded_columns <- unique(c(1, excluded_columns))
  
  # Ignore ID and properties
  raw_data <- raw_data[, -excluded_columns]
  # Compute feature indices, accounting for the offset of ID and properties
  feature_indices <- 1:ncol(raw_data) 
  # Extract features (column names)
  features <- colnames(raw_data)
  raw_data_normalised <- foreach::foreach(i = feature_indices,
                                  .combine = rbind) %dopar% {
                                    # Create and populate entry for current feature
                                    record <- data.frame( 
                                      index = i,
                                      feature = features[i],
                                      values = raw_data[, i],
                                      flag = "Non-normal",
                                      transf = "",
                                      transf.value = NA,
                                      stringsAsFactors = FALSE
                                    )
                                    
                                    # Verify the current feature has at least 3 non-NA rows
                                    if(sum(is.finite(raw_data[, i]), na.rm = TRUE) > 2) {
                                      # Assess normality of feature before transforming it
                                      pvalue <- shapiro.test(raw_data[, i])[[2]]
                                      if(pvalue <= 0.05) { # Data must be transformed
                                        tmp <- MetaPipe::transform_data(data = raw_data[, i], 
                                                              feature = features[i], 
                                                              alpha = 0.05, 
                                                              index = i, 
                                                              transf_vals = transf_vals, 
                                                              plots_prefix = paste0(plots_dir, "/HIST")
                                        )
                                        
                                        if(length(tmp)) {
                                          tmp$index <- i
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
  return(raw_data_normalised)
}

#' Postprocessing for the normality assessmment of the features
#' @param raw_data data frame containing the raw data
#' @param excluded_columns vector containing the indices of the data set properties, excluded columns
#' @param raw_data_normalised data frame containing the normalised raw data, created with \code{\link{assess_normality}}
#' @param out_prefix prefix for output files, plots, etc.
#' @param pareto_scaling booolean flag to indicate whether or not perform a pareto scaling on the normalised data
#'
#' @return Files containing the normal and non-paramatric normalised data are stored on disk
#' @export
#'
#' @examples
#' example_data <- data.frame(ID = c(1,2,3,4,5), 
#'                            P1 = c("one", "two", "three", "four", "five"), 
#'                            F1 = rnorm(5), 
#'                            F2 = rnorm(5))
#' example_data_normalised <- assess_normality(example_data, c(1, 2))
#' assess_normality_postprocessing(example_data, c(1, 2), example_data_normalised)
#' 
#' @seealso \code{\link{assess_normality}} and \code{\link{assess_normality_stats}}
assess_normality_postprocessing <- function(raw_data, 
                                            excluded_columns,
                                            raw_data_normalised,
                                            out_prefix = "metapipe", 
                                            pareto_scaling = FALSE) {
  # Verify the raw_data_normalised object has the right structure
  expected_columns <- c("index", "feature", "values", "flag", "transf", "transf.value")
  if (!all(expected_columns %in% colnames(raw_data_normalised)))
    stop("raw_data_normalised must be the output of the function assess_normality")
  
  # Separate normal and non-parametric entries from the normalised data
  raw_data_normalised_norm <- raw_data_normalised[raw_data_normalised$flag == "Normal", ]
  raw_data_normalised_non_par <- raw_data_normalised[raw_data_normalised$flag == "Non-normal", ]
  
  # Extract feature names for both normal and non-parametric data
  features_non_par <- unique(as.character(raw_data_normalised_non_par$feature))
  features_norm <- unique(as.character(raw_data_normalised_norm$feature))
  features_non_par_len <- length(features_non_par)
  features_norm_len <- length(features_norm)
  
  # Create new objects with the normalised data for normal features and original
  # raw data for non-parametric features
  raw_data_non_par <- NULL
  raw_data_norm <- NULL
  if (features_non_par_len > 0)
    raw_data_non_par <- raw_data[, features_non_par]
  if (features_norm_len > 0 ) {
    raw_data_norm <- data.frame(matrix(vector(), 
                                       nrow(raw_data_normalised_norm) / features_norm_len, 
                                       features_norm_len,
                                       dimnames = list(c(), features_norm)),
                                stringsAsFactors = FALSE)
    
    for (i in 1:features_norm_len) {
      raw_data_norm[i] <- subset(raw_data_normalised_norm, feature == features_norm[i])$values
    }
  }
  
  # Append excluded columns for scaling, if pareto_scaling == TRUE
  if (!is.null(raw_data_norm)) {
    raw_data_norm <- cbind(raw_data[, excluded_columns], 
                           ifelse(pareto_scaling, 
                                  MetaPipe::paretoscale(raw_data_norm), 
                                  raw_data_norm)
    )
  }
  
  if (!is.null(raw_data_non_par)) {
    raw_data_non_par <- cbind(raw_data[, excluded_columns], 
                              ifelse(pareto_scaling, 
                                     MetaPipe::paretoscale(raw_data_non_par), 
                                     raw_data_non_par)
    )
  }
  
  # Generate basic stats from the normalisation process
  raw_data_rows <- nrow(raw_data)
  norm_features_normalised_count <- nrow(raw_data_normalised_norm[raw_data_normalised_norm$transf == "", ]) / raw_data_rows
  norm_features_count <- nrow(raw_data_normalised_norm) / raw_data_rows
  total_features <- nrow(raw_data_normalised)/raw_data_rows
  transformations <- unique(raw_data_normalised[c("transf", "transf.value")])
  transformations <- transformations[-nrow(transformations), ] # Drop blank transformation, NULL transformation
  sorting <- order(transformations$transf, decreasing = TRUE)
  transformations <- transformations[sorting, ]
  
  # Create data frame containing the stats
  normalisation_stats <- data.frame(key = c("total", "norm_features", "norm_features_normalised"),
                                    values = c(total_features, norm_features_count, norm_features_normalised_count),
                                    stringsAsFactors = FALSE)
  if (nrow(transformations) > 0){
    for(i in 1:nrow(transformations)){
      key <- paste0(transformations$transf[i],"\t",transformations$transf.value[i])
      tmp <- subset(raw_data_normalised_norm, raw_data_normalised_norm$transf == transformations$transf[i])
      tmp <- subset(tmp, transf.value == transformations$transf.value[i])
      value <- nrow(tmp) / raw_data_rows
      normalisation_stats <- rbind(normalisation_stats, c(key, value))
    }
  }
  
  # Write to disk new data structures
  write.csv(raw_data_normalised, file = paste0(out_prefix,"_raw_data_normalised_all.csv"), row.names = FALSE)
  write.csv(raw_data_norm, file = paste0(out_prefix,"_raw_data_norm.csv"), row.names = FALSE)
  write.csv(raw_data_non_par, file = paste0(out_prefix,"_raw_data_non_par.csv"), row.names = FALSE)
  write.csv(normalisation_stats, file = paste0(out_prefix,"_normalisation_stats.csv"), row.names = FALSE)
}

#' Statistics for the normality assessmment of the features. 
#' @description This function uses an output file generated by \code{\link{assess_normality_postprocessing}}, 
#' therefore, it cannot be executed without previously running the postprocessing function.
#'
#' @param out_prefix prefix for output files, plots, etc.
#'
#' @export
#'
#' @examples
#' example_data <- data.frame(ID = c(1,2,3,4,5), 
#'                            P1 = c("one", "two", "three", "four", "five"), 
#'                            F1 = rnorm(5), 
#'                            F2 = rnorm(5))
#' example_data_normalised <- assess_normality(example_data, c(1, 2))
#' assess_normality_postprocessing(example_data, c(1, 2), example_data_normalised)
#' assess_normality_stats()
#' 
#' @seealso \code{\link{assess_normality}} and \code{\link{assess_normality_postprocessing}}
assess_normality_stats <- function(out_prefix = "metapipe") {
  stats_filename <- paste0(out_prefix,"_normalisation_stats.csv")
  if (!file.exists(stats_filename))
    stop(paste0("The file ", stats_filename, " was not found"))
  
  # Loading stats for the normality assessment process
  normalisation_stats <- read.csv(stats_filename)
  total_features <- normalisation_stats[1, 2]
  norm_features_count <- normalisation_stats[2, 2]
  norm_features_normalised_count <- normalisation_stats[3, 2]
  transformations <- normalisation_stats[-c(1:3), ]
  
  # Create message for user with the summary
  msg <- paste0("Total features (excluding all NAs features): \t", total_features)
  msg <- paste0(msg, "\nNormal features (without transformation): \t", (norm_features_count - norm_features_normalised_count))
  msg <- paste0(msg, "\nNormal features (transformed): \t\t\t", norm_features_normalised_count)
  msg <- paste0(msg, "\nTotal Normal features: \t\t\t\t", norm_features_count)
  msg <- paste0(msg, "\nNon-parametric features: \t\t\t", (total_features - norm_features_count),"\n")
  
  if (nrow(transformations) > 0) {
    msg <- paste0(msg, "\nTransformations summary:")
    msg <- paste0(msg, "\n\tf(x)\tValue \t# Features")
    for (i in 1:nrow(transformations)) {
      tmp <- strsplit(as.character(transformations[i, 1]), '\t')[[1]]
      msg <- paste0(msg, "\n\t", tmp[1], "\t", tmp[2],"\t", transformations[i, 2])
    }
  }
  msg <- paste0(msg, "\n\n") # Clean output
  message(msg)
}

qtl_scone <- function(x_data, features, CPUS = 1,  ...) {
  cl <- parallel::makeCluster(ceiling(CPUS))
  doParallel::registerDoParallel(cl)
  x_scone <- foreach(i = 2:ncol(x_norm$pheno),
                     .combine = cbind) %dopar% {
                       # Run single scan
                       scone <- qtl::scanone(x_norm, pheno.col = i,  ...) #model = "normal", method = "hk")
                       if(i == 2) {
                         record <- data.frame(
                           chr = scone$chr,
                           pos = scone$pos,
                           lod = scone$lod,
                           row.names = rownames(scone)
                         )
                         colnames(record)[3] <- features[i]
                       }
                       else {
                         record <- data.frame(data = scone$lod)
                         colnames(record) <- features[i]
                       }
                       record
                     }
  parallel::stopCluster(cl) # Stop cluster
  return(x_scone)
}

qtl_preprocessing <- function(genetic_map, out_prefix = "metapipe") {
#   genetic_map <- read.csv("OriginalMap.csv")
#   colnames(genetic_map)[1] <- "ID"
#   genetic_map$ID <- as.character(genetic_map$ID)
#   
  ## Normal features
  raw_data_norm <- read.csv(paste0(out_prefix,"_raw_data_norm.csv"), stringsAsFactors = FALSE)
#   colnames(raw_data_norm)[1] <- "ID"
#   raw_data_norm$GenoID <- with(raw_data_norm,
#                                gsub(" ", "0", paste0(Generation, "_", sprintf("%3s", as.character(ID))))
#                               )
#   raw_data_norm$ID <- raw_data_norm$GenoID
#   raw_data_norm$GenoID <- NULL
#   
  pheno_norm <- dplyr::inner_join(raw_data_norm, genetic_map, by = "ID")[, colnames(raw_data_norm)]
  pheno_norm$Group <- NULL
  pheno_norm$Generation <- NULL
  geno_norm <- rbind(genetic_map[1:2, ],
                     dplyr::inner_join(pheno_norm, genetic_map, by = "ID")[, colnames(genetic_map)]
                    )
#   
#   
  ## Non-parametric features
  raw_data_non_par <- read.csv(paste0(out_prefix, "_raw_data_non_par.csv"), stringsAsFactors = FALSE)
#   colnames(raw_data_non_par)[1] <- "ID"
#   raw_data_non_par$GenoID <- with(raw_data_non_par,
#                                   gsub(" ", "0", paste0(Generation, "_", sprintf("%3s", as.character(ID))))
                                 # )
#   raw_data_non_par$ID <- raw_data_non_par$GenoID
#   raw_data_non_par$GenoID <- NULL
#   
#   pheno_non_par <- dplyr::inner_join(raw_data_non_par, genetic_map, by = "ID")[, colnames(raw_data_non_par)]
#   pheno_non_par$Group <- NULL
#   pheno_non_par$Generation <- NULL
#   non.parametric.gen <- rbind(genetic_map[1:2, ],
#                               dplyr::inner_join(pheno_non_par, genetic_map, by = "ID")[, colnames(genetic_map)]
#                              )
#   
#   # Clean phenotypic data
#   empty_features_non_par <- sapply(pheno_non_par, function(x) all(is.na(x)) || all(is.infinite(x)))
#   empty_features_norm <- sapply(pheno_norm, function(x) all(is.na(x)) || all(is.infinite(x)))
#   #pheno_non_par.ncols <- ncol(pheno_non_par)
#   #pheno_norm.ncols <- ncol(pheno_norm)
#   pheno_non_par[empty_features_non_par] <- NULL
#   pheno_norm[empty_features_norm] <- NULL
#   
#   if(any(empty_features_non_par)) {
#     print(paste0("The following non-parametric features were removed (NAs):"))
#     print(names(empty_features_non_par)[empty_features_non_par])
#   }
#   
#   if(any(empty_features_norm)) {
#     print(paste0("The following normal features were removed (NAs):"))
#     print(names(empty_features_norm)[empty_features_norm])
#   }
#   
#   # Write genotypic and phenotypic dataset
#   ## Normal features
#   write.csv(geno_norm, file = paste0(OUT_PREFIX,".geno_norm.csv"), row.names=FALSE)
#   write.csv(pheno_norm, file = paste0(OUT_PREFIX,".pheno_norm.csv"), row.names=FALSE)
#   ## Non-parametric features
#   write.csv(non.parametric.gen, file = paste0(OUT_PREFIX,".non.parametric.gen.csv"), row.names=FALSE)
#   write.csv(pheno_non_par, file = paste0(OUT_PREFIX,".pheno_non_par.csv"), row.names=FALSE)
}