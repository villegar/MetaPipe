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

check_alpha <- function(alpha) {
  if (alpha < 0)
    stop("alpha must be non-negative")
  else if (alpha > 1)
    stop("alpha must be a numeric value between 0 and 1")
  else if (!is.numeric(alpha))
    stop("alpha must be a numeric value between 0 and 1")
}

log_transformation <- function(shapiro, data, feature, alpha = 0.05,
                               bases = c(2, exp(1), 3, 4, 5, 6, 7, 8, 9, 10)) {
  check_alpha(alpha = alpha)
  
  record <- data.frame(
    index = i,
    feature = features[i],
    values = meansp[, i],
    flag = "Non-normal",
    transf = "",
    transf.value = 0
  )
  
  for (base in bases) {
    pval <- shapiro.test(log(data, base))[[2]]
    if (pval >= alpha) {
      transformed <- log(data, base)
      if (base == exp(1))
        base <- "e"
      xlab <- paste0("$\\log_{", base, "}(", feature, ")$")
      transformation <- paste0("LOG_", base)
      prefix <- paste0("plots/HIST_", (i - 3), "_", transformation)
      compare_hist(data, transformed, feature, prefix, xlab)
      record$values <- transformed
      record$transf <- "log"
      record$transf.value <- base
      return(record)
    }
  }
  return(root_transformation(shapiro, data, feature))
}

power_transformation <- function(shapiro, data, feature) {
  powers <- c(2, exp(1), 3, 4, 5, 6, 7, 8, 9, 10)
  
  record <- data.frame(
    index = i,
    feature = features[i],
    values = meansp[, i],
    flag = "Non-normal",
    transf = "",
    transf.value = 0
  )
  
  for (p in powers) {
    pval <- shapiro.test(data ^ p)[[2]]
    if (pval >= 0.05) {
      transformed <- data ^ p
      if (p == exp(1))
        p <- "e"
      xlab <- paste0("$(", feature, ")^", p, "$")
      transformation <- paste0("POW_", p)
      prefix <- paste0("plots/HIST_", (i - 3), "_", transformation)
      compare_hist(data, transformed, feature, prefix, xlab)
      record$values <- transformed
      record$transf <- "power"
      record$transf.value <- p
      return(record)
    }
  }
  return(data.frame())
}

root_transformation <- function(shapiro, data, feature) {
  roots <- c(2, exp(1), 3, 4, 5, 6, 7, 8, 9, 10)
  
  record <- data.frame(
    index = i,
    feature = features[i],
    values = meansp[, i],
    flag = "Non-normal",
    transf = "",
    transf.value = 0
  )
  
  for (r in roots) {
    pval <- shapiro.test(data ^ (1 / r))[[2]]
    if (pval >= 0.05) {
      transformed <- data ^ (1 / r)
      if (r == exp(1))
        r <- "e"
      xlab <- paste0("$\\sqrt[", r, "]{", feature, "}$")
      transformation <- paste0("ROOT_", r)
      prefix <- paste0("plots/HIST_", (i - 3), "_", transformation)
      compare_hist(data, transformed, feature, prefix, xlab)
      record$values <- transformed
      record$transf <- "root"
      record$transf.value <- r
      return(record)
    }
  }
  return(power_transformation(shapiro, data, feature))
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
