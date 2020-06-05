#' Function to save a graphical object to disk in 
#' Portable Network Graphics (PNG) format.
#'
#' @param myPlot graphical object
#' @param name output name with or without path
#' @param width width in inches
#' @param height height in inches
#'
#' @export
#'
#' @examples
#' savePlot(hist(rnorm(100), main = "Histogram of Normal Distribution"), 
#' "hist")
#' 
#' @seealso \code{\link{savePlotPDF}} and \code{\link{savePlotTIFF}}
savePlot <- function(myPlot,name, width = 6, height = 6) {
  grDevices::png(paste0(name,".png"), width = width, height = height, units = 'in', res = 300, type = "cairo")
  print(myPlot)
  grDevices::dev.off()
}

#' Function to save a graphical object to disk in 
#' Tagged Image File Format (TIFF) format.
#'
#' @param myPlot graphical object
#' @param name output name with or without path
#' @param width width in inches
#' @param height height in inches
#'
#' @export
#'
#' @examples
#' savePlotTIFF(hist(rnorm(100), main = "Histogram of Normal Distribution"), 
#' "hist")
#' 
#' @seealso \code{\link{savePlotPDF}} and \code{\link{savePlot}}
savePlotTIFF <- function(myPlot,name, width = 6, height = 6) {
  grDevices::tiff(paste0(name,".tiff"), width = width, height = height, units = 'in', res = 300, type = "cairo")
  print(myPlot)
  grDevices::dev.off()
}

#' Function to save a graphical object to disk in 
#' Portable Document Format (PDF) format.
#'
#' @param myPlot graphical object
#' @param name output name with or without path
#' @param width width in inches
#' @param height height in inches
#'
#' @export
#'
#' @examples
#' savePlotPDF(hist(rnorm(100), main = "Histogram of Normal Distribution"), 
#' "hist")
#' 
#' @seealso \code{\link{savePlot}} and \code{\link{savePlotTIFF}}
savePlotPDF <- function(myPlot,name, width = 6, height = 6) {
  grDevices::pdf(paste0(name,".pdf"), width = width, height = height)
  print(myPlot)
  grDevices::dev.off()
}

ggplot_save <- function(myPlot,name){
  R.devices::suppressGraphics({
    ggplot2::ggsave(
      paste0(name,".png"),
      plot   = myPlot,
      device = 'png',
      width  = 6,
      height = 6,
      dpi    = 300,
      type   = "cairo"
    )
  })
}

# This function generates a two histograms to compare two datasets
## original = original data
## transformed = transformed data
## feature = data feature
## name.prefix = prefix for the file name
## xlab = label for x-axis
#' Title
#'
#' @param original Original data
#' @param transformed Transformed data
#' @param feature Feature name
#' @param name.prefix File prefix
#' @param xlab x-axis label
#'
#' @export
#'
#' @examples
#' norm_dist <- rnorm(100)
#' norm_dist_transformed <- norm_dist^2
#' compare_hist(norm_dist, norm_dist_transformed, "XYZ", "xyz_hist","")
compare_hist <- function(original,transformed,feature,name.prefix,xlab){
  ALPHA <- 1
  BINS <- 20
  histogram <- data.frame(
    original = original, 
    transformed = transformed
  )
  original.plot <- ggplot2::ggplot(data = histogram, ggplot2::aes(original)) +
    ggplot2::geom_histogram(alpha = ALPHA, ggplot2::aes(y = ..count..), position = 'identity', bins = BINS, col = "black", fill = "#FFDF01") + #"#B2182B") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1)) + 
    ggplot2::labs(title=paste("Feature",feature), x='', y='') +
    ggplot2::xlab(paste0(feature))
  
  transformed.plot <- ggplot2::ggplot(data = histogram, ggplot2::aes(transformed)) +
    ggplot2::geom_histogram(alpha = ALPHA, ggplot2::aes(y = ..count..), position = 'identity', bins = BINS, col = "black", fill = "#0057A7") + #"#2166AC") +
    # geom_density(ggplot2::aes(y=..density..)) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1)) + 
    ggplot2::labs(x='', y='') +
    ggplot2::xlab(latex2exp::TeX(xlab))
  
  ggplot_save(grid::grid.draw(rbind(ggplot2::ggplotGrob(original.plot),ggplot2::ggplotGrob(transformed.plot), size = "last")),
              paste0(name.prefix,"_",feature))
}

# This function generates a single histogram
## data = input data
## feature = data feature
## name.prefix = prefix for the file name
## xlab = label for x-axis
generate.histogram <- function(data,feature,name.prefix,xlab){
  ALPHA <- 1
  BINS <- 20
  histogram <- data.frame(original = data)
  myPlot <- ggplot2::ggplot(data = histogram, ggplot2::aes(histogram$original)) +
    ggplot2::geom_histogram(alpha = ALPHA, ggplot2::aes(y = ggplot2::stat_count()), position = 'identity', bins = BINS, col = "black", fill = "#7FCDBB") + #"#90C978") +
    # geom_density(ggplot2::aes(y=..density..)) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1)) + 
    ggplot2::labs(title=paste("Feature",feature), x='', y='') +
    ggplot2::xlab(paste0(feature))
  ggplot_save(myPlot,paste0(name.prefix,"_",feature))
}
