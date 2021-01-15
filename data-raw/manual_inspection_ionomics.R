
mad_outliers <- function(data, mult = 1) {
  # lower_bound <- median(data) - 3 * mad(data)
  # upper_bound <- median(data) + 3 * mad(data)
  # idx <- which(data < mult * lower_bound | data > mult * upper_bound)
  # po <- boxplot.stats(data)$out
  # lower_bound <- median(po) - 3 * mad(po)
  # upper_bound <- median(po) + 3 * mad(po)
  # idx <- which(po < mult * lower_bound | po > mult * upper_bound)
  idx <- which(data %in% boxplot.stats(data)$out)
  if (length(idx) > 0) {
    col <- rep("black", length(data))
    col[idx] <- "red"
    plot(data, col = col)
    df <- data.frame(pos = idx, val = data[idx])
    print(df)
    # message(paste0("IQR: ", IQR(data, na.rm = TRUE)))
    # message(paste0(idx, collapse = ", "))
    # print(data[idx])
    return(TRUE)
  }
  return(FALSE)
}
# sapply(ionomics[, -c(1:2)], mad_outliers)
mad_outliers(ionomics[, 4])
mad_outliers(ionomics[, 5])
mad_outliers(ionomics[, 6])
ionomics[ionomics[, 6] > 0.10, 6] <- NA
mad_outliers(ionomics[, 7])
mad_outliers(ionomics[, 8])
ionomics[ionomics[, 8] > 50, 8] <- NA # Significant
mad_outliers(ionomics[, 9])
mad_outliers(ionomics[, 10])
ionomics[ionomics[, 10] > 40, 10] <- NA
mad_outliers(ionomics[, 11])
mad_outliers(ionomics[, 12])
ionomics[ionomics[, 12] > 100, 12] <- NA
mad_outliers(ionomics[, 13])
mad_outliers(ionomics[, 14])
ionomics[ionomics[, 14] > 0.4, 14] <- NA
mad_outliers(ionomics[, 15])
ionomics[ionomics[, 15] > 1.5, 15] <- NA
mad_outliers(ionomics[, 16])
mad_outliers(ionomics[, 17])
mad_outliers(ionomics[, 18])
mad_outliers(ionomics[, 19])
ionomics[abs(ionomics[19]) > 1.5, 19] <- NA
mad_outliers(ionomics[, 20])
mad_outliers(ionomics[, 21])