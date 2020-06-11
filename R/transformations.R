# Log Transformation and Pareto scaling
## Pareto scaling function
#This function is adapted from Stephen C. Grace and Dane A. Hudson 
paretoscale <- function(z) {
  colmean <- apply(z, 2, mean)
  colsd <- apply(z, 2, sd)
  colsqrtsd <- sqrt(colsd)
  rv <- sweep(z, 2, colmean, "-")
  rv <- sweep(rv, 2, colsqrtsd, "/")
  return(rv)
}

#' Check significance level (alpha) validity 
#'
#' @param alpha significance level
#'
#' @return Nothing if alpha is valid, otherwise stop execution
#' @export
#'
#' @examples
#' check_alpha(0.5)
#' \dontrun{
#' check_alpha(-0.5)
#' check_alpha(1.5)
#' check_alpha("1.5")
#' }
check_alpha <- function(alpha) {
  if (!is.numeric(alpha))
    stop("alpha must be a numeric value")
  else if (alpha < 0)
    stop("alpha must be non-negative")
  else if (alpha > 1)
    stop("alpha must be a numeric value between 0 and 1")
}

check_transformation <- function(ref, new, transf) {
  if (new < ref) { # Compare a new p-value with a reference (original).
    warning(paste0(transf, " transformation does not normalise the data."))
    return(TRUE)
  }
  return(FALSE)
}

#' Normalise data with a log transformation
#'
#' @param data original data
#' @param feature feature name
#' @param alpha significance level
#' @param transf transformation values
#' @param plots_prefix prefix for plots with or without path
#' @param digits significant digits to compare p-values of transformations
#'
#' @return data structure containing the normalised data, NULL if no 
#' transformation was performed
#' @export
#'
#' @examples
#' set.seed(123)
#' data <- rnorm(100)
#' log_transformation(2 ^ data, "EXP_2")
log_transformation <- function(data, feature, alpha = 0.05,
                               transf = c(2, exp(1), 3, 4, 5, 6, 7, 8, 9, 10),
                               plots_prefix = "HIST",
                               digits = 6) {
  check_alpha(alpha = alpha) # Check the significance level
  
  ref_pval <- shapiro.test(data)[[2]] # Obtain reference p-value
  
  # Data frame to hold the outcome of the transformation
  record <- data.frame(
    flag = "Non-normal",
    transf = "",
    transf.value = 0
  )
  
  pvals <- data.frame(matrix(vector(), 1, length(transf)))
  for (k in 1:ncol(pvals)) {
    suppressWarnings({
      pvals[1, k] <- shapiro.test(log(data, transf[k]))[[2]]
    })
  }
  
  # Obtain transformation with largest p-value
  max_pval_idx <- unname(which.max(round(pvals, digits)))
  max_pval <- max(pvals, na.rm = TRUE)
  
  # Verify if a transformation normalised the data
  if (check_transformation(alpha, max_pval, "Log"))
    return(NULL) 
  if (check_transformation(ref_pval, max_pval, "Log"))
    return(NULL) 
  
  base <- transf[max_pval_idx]
  transformed <- log(data, base)
  
  base <- ifelse(base == exp(1), "e", base)
  
  xlab <- paste0("$\\log_{", base, "}(", feature, ")$")
  transformation <- paste0("LOG_", base)
  prefix <- paste0(plots_prefix, "_", transformation)
  compare_hist(data, transformed, feature, prefix, xlab)
  record$flag = "Normal"
  record$transf <- "log"
  record$transf.value <- base
  return(record)
}

#' Normalise data with a power transformation
#'
#' @param data original data
#' @param feature feature name
#' @param alpha significance level
#' @param transf transformation values
#' @param plots_prefix prefix for plots with or without path
#' @param digits significant digits to compare p-values of transformations
#'
#' @return data structure containing the normalised data, NULL if no 
#' transformation was performed
#' @export
#'
#' @examples
#' set.seed(123)
#' data <- rnorm(100, 5)
#' power_transformation(sqrt(data), "ROOT_2")
power_transformation <- function(data, feature, alpha = 0.05,
                                 transf = c(2, exp(1), 3, 4, 5, 6, 7, 8, 9, 10),
                                 plots_prefix = "HIST",
                                 digits = 6) {
  check_alpha(alpha = alpha) # Check the significance level
  
  ref_pval <- shapiro.test(data)[[2]] # Obtain reference p-value
  
  # Data frame to hold the outcome of the transformation
  record <- data.frame(
    flag = "Non-normal",
    transf = "",
    transf.value = 0
  )
  
  pvals <- data.frame(matrix(vector(), 1, length(transf)))
  for (k in 1:ncol(pvals)) {
    suppressWarnings({
      pvals[1, k] <- shapiro.test(data ^ transf[k])[[2]]
    })
  }
  
  # Obtain transformation with largest p-value
  max_pval_idx <- unname(which.max(round(pvals, digits)))
  max_pval <- max(pvals, na.rm = TRUE)
  
  # Verify if a transformation normalised the data
  if (check_transformation(alpha, max_pval, "Power"))
    return(NULL) 
  if (check_transformation(ref_pval, max_pval, "Power"))
    return(NULL) 
  
  power <- transf[max_pval_idx]
  transformed <- data ^ power
  
  power <- ifelse(power == exp(1), "e", power)
  
  xlab <- paste0("$(", feature, ")^", power, "$")
  transformation <- paste0("POW_", power)
  prefix <- paste0(plots_prefix, "_", transformation)
  compare_hist(data, transformed, feature, prefix, xlab)
  record$flag = "Normal"
  record$transf <- "power"
  record$transf.value <- power
  return(record)
}

#' Normalise data with a root transformation
#'
#' @param data original data
#' @param feature feature name
#' @param alpha significance level
#' @param transf transformation values
#' @param plots_prefix prefix for plots with or without path
#' @param digits significant digits to compare p-values of transformations
#'
#' @return data structure containing the normalised data, NULL if no 
#' transformation was performed
#' @export
#'
#' @examples
#' set.seed(123)
#' data <- rnorm(100, 5)
#' root_transformation(data ^ 2, "EXP_2")
root_transformation <- function(data, feature, alpha = 0.05,
                                transf = c(2, exp(1), 3, 4, 5, 6, 7, 8, 9, 10),
                                plots_prefix = "HIST",
                                digits = 6) {
  
  check_alpha(alpha = alpha) # Check the significance level
  
  ref_pval <- shapiro.test(data)[[2]] # Obtain reference p-value
  
  # Data frame to hold the outcome of the transformation
  record <- data.frame(
    flag = "Non-normal",
    transf = "",
    transf.value = 0
  )
  
  pvals <- data.frame(matrix(vector(), 1, length(transf)))
  for (k in 1:ncol(pvals)) {
    suppressWarnings({
      pvals[1, k] <- shapiro.test(data ^ (1 / transf[k]))[[2]]
    })
  }
  
  # Obtain transformation with largest p-value
  max_pval_idx <- unname(which.max(round(pvals, digits)))
  max_pval <- max(pvals, na.rm = TRUE)
  
  # Verify if a transformation normalised the data
  if (check_transformation(alpha, max_pval, "Root"))
    return(NULL) 
  if (check_transformation(ref_pval, max_pval, "Root"))
    return(NULL) 
  
  root <- transf[max_pval_idx]
  transformed <- data ^ (1 / root)
  
  root <- ifelse(root == exp(1), "e", root)
  
  xlab <- paste0("$\\sqrt[", root, "]{", feature, "}$")
  transformation <- paste0("ROOT_", root)
  prefix <- paste0(plots_prefix, "_", transformation)
  compare_hist(data, transformed, feature, prefix, xlab)
  record$flag = "Normal"
  record$transf <- "root"
  record$transf.value <- root
  return(record)
}

transform_data <- function(shapiro,
                           data,
                           feature,
                           index,
                           offset = 3,
                           plots.directory = "plots",
                           transformation.values = 
                             c(2, exp(1), 3, 4, 5, 6, 7, 8, 9, 10)) {
  #bases <- c(2,exp(1),3,4,5,6,7,8,9,10)
  #powers <- c(2,exp(1),3,4,5,6,7,8,9,10)
  #roots <- c(2,exp(1),3,4,5,6,7,8,9,10)
  
  record <- data.frame(
    index = index,
    feature = feature,
    values = data,
    flag = "Non-normal",
    transf = "",
    transf.value = 0
  )
  
  pvals <- data.frame(matrix(vector(), 3, length(transformation.values)))
  for (k in 1:ncol(pvals)) {
    suppressWarnings({
      pvals[1, k] <- shapiro.test(log(data, transformation.values[k]))[[2]]
      pvals[2, k] <- shapiro.test(data ^ transformation.values[k])[[2]]
      pvals[3, k] <- shapiro.test(data ^ (1 / transformation.values[k]))[[2]]
    })
  }
  max.pval <- max(pvals, na.rm = TRUE)
  max.pval.index <- which(pvals == max.pval, arr.ind = TRUE)
  
  transf <- max.pval.index[1]
  transf.value.index <- max.pval.index[2]
  
  if (max.pval < 0.05) # Verify if a transformation normalized the data
    return(data.frame())
  
  if (transf == 1) { # Log transformation
    base <- transformation.values[transf.value.index]
    transformed <- log(data, base)
    if (base == exp(1))
      base <- "e"
    xlab <- paste0("$\\log_{", base, "}(", feature, ")$")
    transformation <- paste0("LOG_", base)
    prefix <-
      paste0(plots.directory,
             "/HIST_",
             (index - offset),
             "_",
             transformation)
    compare_hist(data, transformed, feature, prefix, xlab)
    record$values <- transformed
    record$transf <- "log"
    record$transf.value <- base
    return(record)
  }
  else if (transf == 2) { # Power transformation
    p <- transformation.values[transf.value.index]
    transformed <- data ^ p
    if (p == exp(1))
      p <- "e"
    xlab <- paste0("$(", feature, ")^", p, "$")
    transformation <- paste0("POW_", p)
    prefix <-
      paste0(plots.directory,
             "/HIST_",
             (index - offset),
             "_",
             transformation)
    compare_hist(data, transformed, feature, prefix, xlab)
    record$values <- transformed
    record$transf <- "power"
    record$transf.value <- p
    return(record)
  }
  else { # Root transformation
    r <- transformation.values[transf.value.index]
    transformed <- data ^ (1 / r)
    if (r == exp(1))
      r <- "e"
    xlab <- paste0("$\\sqrt[", r, "]{", feature, "}$")
    transformation <- paste0("ROOT_", r)
    prefix <-
      paste0(plots.directory,
             "/HIST_",
             (index - offset),
             "_",
             transformation)
    compare_hist(data, transformed, feature, prefix, xlab)
    record$values <- transformed
    record$transf <- "root"
    record$transf.value <- r
    return(record)
  }
}
