#' Pareto scaling
#' 
#' Perform a pareto scaling.
#' 
#' @importFrom stats sd
#' @param z Original 2-dimensional data.
#'
#' @return Scaled data.
#'
#' @examples
#' set.seed(123)
#' data <- matrix(rnorm(100, 5), ncol = 2)
#' data_new <- MetaPipe:::paretoscale(data)
#' 
#' @keywords internal
#' @noRd
#' 
#' @references 
#' Grace, Stephen & Hudson, Dane. (2016). Processing and Visualization of 
#' Metabolomics Data Using R. \doi{10.5772/65405}. 
paretoscale <- function(z) {
  colmean <- apply(z, 2, mean)
  colsd <- apply(z, 2, sd)
  colsqrtsd <- sqrt(colsd)
  rv <- sweep(z, 2, colmean, "-")
  rv <- sweep(rv, 2, colsqrtsd, "/")
  return(rv)
}

#' Check significance level
#' 
#' Check significance level (\code{alpha}) validity.
#'
#' @param alpha Significance level.
#'
#' @return Nothing if alpha is valid, otherwise stop execution.
#'
#' @examples
#' MetaPipe:::check_alpha(0.5)
#' \dontrun{
#'     MetaPipe:::check_alpha(-0.5)
#'     MetaPipe:::check_alpha(1.5)
#'     MetaPipe:::check_alpha("1.5")
#' }
#' 
#' @keywords internal
#' @noRd
check_alpha <- function(alpha) {
  if (!(typeof(alpha) %in% c("integer", "double")))
    stop("alpha must be a numeric value")
  else if (alpha < 0)
    stop("alpha must be non-negative")
  else if (alpha > 1)
    stop("alpha must be a numeric value between 0 and 1")
}

#' Check transformation
#' 
#' Check whether or not the transformation normalises the data.
#'
#' @param ref Reference p-value.
#' @param new New p-value.
#' @param transf Transformation name (e.g. Log).
#' @param msg Message template.
#'
#' @return \code{TRUE} if the transformation normalised the data, \code{FALSE} 
#' otherwise.
#'
#' @examples
#' MetaPipe:::check_transformation(0.05, 0.045, "Log")
#' MetaPipe:::check_transformation(0.05, 0.1, "Log")
#' 
#' @keywords internal
#' @noRd
check_transformation <- 
  function(ref, 
           new, 
           transf, 
           msg = " transformation does not normalise the data.") {
  if (new < ref) { # Compare a new p-value with a reference (original).
    message(paste0("\n", transf, msg))
    return(TRUE)
  }
  return(FALSE)
}

#' Logarithmic transformation
#' 
#' Normalise data with a log transformation.
#' 
#' @param data Original data.
#' @param trait Trait name.
#' @param alpha Significance level.
#' @param transf_vals Transformation values.
#' @param plots_prefix Prefix for plots with or without path.
#' @param digits Significant digits to compare p-values of transformations.
#' @param plot Boolean flag to indicate whether or not to plot a histogram 
#'     comparing the original and transformed/normalised data.
#' @param quiet Boolean flag to hide status messages. 
#'
#' @return Data structure containing the normalised data, original data if no 
#' transformation was performed.
#' 
#' @export
#'
#' @examples
#' set.seed(123)
#' data <- rnorm(100, 5)
#' transformed <- MetaPipe::log_transformation(2 ^ data, "EXP_2")
#' 
#' # Clean up example outputs
#' MetaPipe:::tidy_up("HIST")
#' 
#' @family Transformation functions
log_transformation <- function(data, 
                               trait = "DATA", 
                               alpha = 0.05,
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
                               plots_prefix = "HIST",
                               digits = 6, 
                               plot = TRUE,
                               quiet = FALSE) {
  check_alpha(alpha = alpha) # Check the significance level
  
  ref_pval <- shapiro.test(data)[[2]] # Obtain reference p-value
  
  # # Data frame to hold the outcome of the transformation
  # record <- data.frame(
  #   flag = "Skewed",
  #   transf = "",
  #   transf_val = 0
  # )
  
  pvals <- data.frame(matrix(vector(), 1, length(transf_vals)))
  for (k in 1:ncol(pvals)) {
    suppressWarnings({
      pvals[1, k] <- shapiro.test(log(data, transf_vals[k]))[[2]]
    })
  }
  
  # Obtain transformation with largest p-value
  max_pval_idx <- unname(which.max(round(pvals, digits)))
  max_pval <- max(pvals, na.rm = TRUE)
  
  # Verify if a transformation normalised the data
  if (check_transformation(alpha, max_pval, "Log") ||
      check_transformation(ref_pval, max_pval, "Log")) {
    warning("No logarithmic transformation was found to transform/normalise ",
            "the data.")
    return(data)
  }
  
  base <- transf_vals[max_pval_idx]
  transformed <- log(data, base)
  
  base <- ifelse(base == exp(1), "e", base)
  
  if (plot) {
    xlab <- paste0("$\\log_{", base, "}(", trait, ")$")
    transformation <- paste0("LOG_", base)
    prefix <- paste0(plots_prefix, "_", transformation)
    compare_hist(data, transformed, trait, prefix, xlab)
  }
  # record$flag = "Normal"
  # record$transf <- "log"
  # record$transf_val <- base
  if (!quiet)
    message("\nThe data was normalised by a logarithm of base ", base)
  return(transformed)
}

#' Power transformation
#' 
#' Normalise data with a power transformation.
#'
#' @param data Original data.
#' @param trait Trait name.
#' @param alpha Significance level.
#' @param transf_vals Transformation values.
#' @param plots_prefix Prefix for plots with or without path.
#' @param digits Significant digits to compare p-values of transformations.
#' @param plot Boolean flag to indicate whether or not to plot a histogram 
#'     comparing the original and transformed/normalised data.
#' @param quiet Boolean flag to hide status messages. 
#'
#' @return Data structure containing the normalised data, original data if no 
#' transformation was performed.
#' 
#' @export
#'
#' @examples
#' set.seed(123)
#' data <- rnorm(100, 5)
#' transformed <- MetaPipe::power_transformation(sqrt(data), "ROOT_2")
#' 
#' # Clean up example outputs
#' MetaPipe:::tidy_up("HIST")
#' 
#' @family Transformation functions
power_transformation <- function(data, 
                                 trait = "DATA", 
                                 alpha = 0.05,
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
                                 plots_prefix = "HIST",
                                 digits = 6, 
                                 plot = TRUE,
                                 quiet = FALSE) {
  check_alpha(alpha = alpha) # Check the significance level
  
  ref_pval <- shapiro.test(data)[[2]] # Obtain reference p-value
  
  # # Data frame to hold the outcome of the transformation
  # record <- data.frame(
  #   flag = "Skewed",
  #   transf = "",
  #   transf_val = 0
  # )
  
  pvals <- data.frame(matrix(vector(), 1, length(transf_vals)))
  for (k in 1:ncol(pvals)) {
    suppressWarnings({
      pvals[1, k] <- shapiro.test(data ^ transf_vals[k])[[2]]
    })
  }
  
  # Obtain transformation with largest p-value
  max_pval_idx <- unname(which.max(round(pvals, digits)))
  max_pval <- max(pvals, na.rm = TRUE)
  
  # Verify if a transformation normalised the data
  if (check_transformation(alpha, max_pval, "Power") || 
      check_transformation(ref_pval, max_pval, "Power")){
    warning("No power transformation was found to transform/normalise ",
            "the data.")
    return(data)
  }
  
  power <- transf_vals[max_pval_idx]
  transformed <- data ^ power
  
  power <- ifelse(power == exp(1), "e", power)
  
  if (plot) {
    xlab <- paste0("$(", trait, ")^", power, "$")
    transformation <- paste0("POW_", power)
    prefix <- paste0(plots_prefix, "_", transformation)
    compare_hist(data, transformed, trait, prefix, xlab)
  }
  # record$flag = "Normal"
  # record$transf <- "power"
  # record$transf_val <- power
  if (!quiet)
    message("\nThe data was normalised by a power of ", power)
  return(transformed)
}

#' Root transformation
#' 
#' Normalise data with a root transformation.
#'
#' @param data Original data.
#' @param trait Trait name.
#' @param alpha Significance level.
#' @param transf_vals Transformation values.
#' @param plots_prefix Prefix for plots with or without path.
#' @param digits Significant digits to compare p-values of transformations.
#' @param plot Boolean flag to indicate whether or not to plot a histogram 
#'     comparing the original and transformed/normalised data.
#' @param quiet Boolean flag to hide status messages. 
#'
#' @return Data structure containing the normalised data, original data if no 
#' transformation was performed.
#' @export
#'
#' @examples
#' set.seed(123)
#' data <- rnorm(100, 5)
#' transformed <- MetaPipe::root_transformation(data ^ 2, "EXP_2")
#' 
#' # Clean up example outputs
#' MetaPipe:::tidy_up("HIST")
#' 
#' @family Transformation functions
root_transformation <- function(data, 
                                trait = "DATA", 
                                alpha = 0.05,
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
                                plots_prefix = "HIST",
                                digits = 6, 
                                plot = TRUE,
                                quiet = FALSE) {
  
  check_alpha(alpha = alpha) # Check the significance level
  
  ref_pval <- shapiro.test(data)[[2]] # Obtain reference p-value
  
  # # Data frame to hold the outcome of the transformation
  # record <- data.frame(
  #   flag = "Skewed",
  #   transf = "",
  #   transf_val = 0
  # )
  
  pvals <- data.frame(matrix(vector(), 1, length(transf_vals)))
  for (k in 1:ncol(pvals)) {
    suppressWarnings({
      pvals[1, k] <- shapiro.test(data ^ (1 / transf_vals[k]))[[2]]
    })
  }
  
  # Obtain transformation with largest p-value
  max_pval_idx <- unname(which.max(round(pvals, digits)))
  max_pval <- max(pvals, na.rm = TRUE)
  
  # Verify if a transformation normalised the data
  if (check_transformation(alpha, max_pval, "Root") ||
      check_transformation(ref_pval, max_pval, "Root")){
    warning("No root transformation was found to transform/normalise ",
            "the data.")
    return(data)
  }
  
  root <- transf_vals[max_pval_idx]
  transformed <- data ^ (1 / root)
  
  root <- ifelse(root == exp(1), "e", root)
  
  if (plot) {
    xlab <- paste0("$\\sqrt[", root, "]{", trait, "}$")
    transformation <- paste0("ROOT_", root)
    prefix <- paste0(plots_prefix, "_", transformation)
    compare_hist(data, transformed, trait, prefix, xlab)
  }
  # record$flag = "Normal"
  # record$transf <- "root"
  # record$transf_val <- root
  if (!quiet)
    message("\nThe data was normalised by a root of ", root)
  return(transformed)
}

#' Normalise data
#' 
#' Normalise data using different methods: log, power, and root, with a 
#' different number of parameters to find out which one transforms the data
#' into a normal-ish set. 
#'
#' @param data Original data.
#' @param trait Trait name.
#' @param alpha Significance level.
#' @param index Index of the current trait.
#' @param transf_vals Transformation values.
#' @param plots_prefix Prefix for plots with or without path.
#' @param digits Significant digits to compare p-values of transformations.
#'
#' @return Transformed record.
#' @export
#'
#' @examples
#' set.seed(123)
#' data <- rnorm(100, 5)
#' out_exp2 <- MetaPipe::transform_data(2 ^ data, "EXP_2")
#' out_root2 <- MetaPipe::transform_data(sqrt(data), "ROOT_2")
#' out_pow2 <- MetaPipe::transform_data(data ^ 2, "POW_2")
#' 
#' knitr::kable(head(out_exp2))
#' knitr::kable(head(out_root2))
#' knitr::kable(head(out_pow2))
#' 
#' # Clean up example outputs
#' MetaPipe:::tidy_up("HIST")
transform_data <- function(data,
                           trait = "DATA",
                           alpha = 0.05,
                           index = NULL,
                           transf_vals = c(2, exp(1), 3, 4, 5, 6, 7, 8, 9, 10),
                           plots_prefix = "HIST",
                           digits = 6) {
  
  check_alpha(alpha = alpha) # Check the significance level
  
  ref_pval <- shapiro.test(data)[[2]] # Obtain reference p-value
  
  trait <- ifelse(is.null(trait), "", trait)
  index <- ifelse(is.null(index), "", paste0(index, "_"))
  
  record <- data.frame(
    index = index,
    trait = trait,
    values = data,
    flag = "Skewed",
    transf = "",
    transf_val = 0,
    stringsAsFactors = FALSE
  )
  
  pvals <- data.frame(matrix(vector(), 3, length(transf_vals)))
  for (k in seq_len(ncol(pvals))) {
    suppressWarnings({
      pvals[1, k] <- tryCatch(
        round(shapiro.test(log(data, transf_vals[k]))[[2]], digits),
        error = function(e) NA) # capture conversion errors
      pvals[2, k] <- tryCatch(
        round(shapiro.test(data ^ transf_vals[k])[[2]], digits),
        error = function(e) NA)  # capture conversion errors
      pvals[3, k] <- tryCatch(
        round(shapiro.test(data ^ (1 / transf_vals[k]))[[2]], digits),
        error = function(e) NA) # capture conversion errors
    })
  }
  
  # Obtain transformation with largest p-value
  # max_pval_idx <- unname(which.max(round(pvals, digits)))
  # max_pval <- max(pvals, na.rm = TRUE)
  max_pval <- max(pvals, na.rm = TRUE)
  max_pval_idx <- which(pvals == max_pval, arr.ind = TRUE)
  
  transf <- max_pval_idx[1]
  transf_val_idx <- max_pval_idx[2]
  
  # Verify if a transformation normalised the data
  if (check_transformation(alpha, 
                           max_pval, 
                           "", 
                           "No transformation normalised the data.") || 
      check_transformation(ref_pval, 
                           max_pval, 
                           "", 
                           "No transformation normalised the data."))
    return(NULL)
  
  record$flag <- "Normal"
  if (transf == 1) { # Log transformation
    base <- transf_vals[transf_val_idx]
    transformed <- log(data, base)
    base <- ifelse(base == exp(1), "e", base)

    xlab <- paste0("$\\log_{", base, "}(", trait, ")$")
    transformation <- paste0("LOG_", base)
    prefix <-
      paste0(plots_prefix,
             "_",
             index,
             transformation)
    compare_hist(data, transformed, trait, prefix, xlab)
    record$values <- transformed
    record$transf <- "log"
    record$transf_val <- base
    return(record)
  } else if (transf == 2) { # Power transformation
    power <- transf_vals[transf_val_idx]
    transformed <- data ^ power
    power <- ifelse(power == exp(1), "e", power)

    xlab <- paste0("$(", trait, ")^", power, "$")
    transformation <- paste0("POW_", power)
    prefix <-
      paste0(plots_prefix,
             "_",
             index,
             transformation)
    compare_hist(data, transformed, trait, prefix, xlab)
    record$values <- transformed
    record$transf <- "power"
    record$transf_val <- power
    return(record)
  } else { # Root transformation
    root <- transf_vals[transf_val_idx]
    transformed <- data ^ (1 / root)
    root <- ifelse(root == exp(1), "e", root)
    
    xlab <- paste0("$\\sqrt[", root, "]{", trait, "}$")
    transformation <- paste0("ROOT_", root)
    prefix <-
      paste0(plots_prefix,
             "_",
             index,
             transformation)
    compare_hist(data, transformed, trait, prefix, xlab)
    record$values <- transformed
    record$transf <- "root"
    record$transf_val <- root
    return(record)
  }
}
