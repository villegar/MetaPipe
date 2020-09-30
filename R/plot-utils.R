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
#' \dontrun{
#'     save_plot(hist(rnorm(100), 
#'               main = "Histogram of Normal Distribution"), 
#'               "hist")
#' }
#' 
#' @seealso \code{\link{save_plotPDF}} and \code{\link{save_plotTIFF}}
save_plot <- function(myPlot, name, width = 6, height = 6) {
  grDevices::png(paste0(name, ".png"), 
                 width = width, 
                 height = height, 
                 units = "in", 
                 res = 300, 
                 type = "cairo")
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
#' \dontrun{
#'     save_plotTIFF(hist(rnorm(100), 
#'                   main = "Histogram of Normal Distribution"), 
#'                   "hist")
#' }
#' 
#' @seealso \code{\link{save_plotPDF}} and \code{\link{save_plot}}
save_plotTIFF <- function(myPlot, name, width = 6, height = 6) {
  grDevices::tiff(paste0(name, ".tiff"), 
                  width = width, 
                  height = height, 
                  units = "in", 
                  res = 300, 
                  type = "cairo")
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
#' \dontrun{
#'     save_plotPDF(hist(rnorm(100), 
#'                  main = "Histogram of Normal Distribution"), 
#'                  "hist")
#' }
#' 
#' @seealso \code{\link{save_plot}} and \code{\link{save_plotTIFF}}
save_plotPDF <- function(myPlot, name, width = 6, height = 6) {
  grDevices::pdf(paste0(name, ".pdf"),
                 width = width,
                 height = height)
  print(myPlot)
  grDevices::dev.off()
}

#' Function to save a graphical object generated with \code{ggplot2} to disk in 
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
#' \dontrun{
#'     myplot <- ggplot2::qplot(rnorm(100))
#'     ggplot_save(myplot, "hist")
#' }
ggplot_save <- function(myPlot, name, width = 6, height = 6){
  R.devices::suppressGraphics({
    ggplot2::ggsave(
      paste0(name, ".png"),
      plot   = myPlot,
      device = "png",
      width  = width,
      height = height,
      dpi    = 300,
      type   = "cairo"
    )
  })
}

#' This function generates a two histograms to compare two datasets
#'
#' @param original Original data
#' @param transformed Transformed data
#' @param feature Feature name
#' @param prefix File prefix
#' @param xlab x-axis label
#'
#' @export
#'
#' @examples
#' \dontrun{
#'     norm_dist <- rnorm(100)
#'     norm_dist_transformed <- norm_dist^2
#'     compare_hist(norm_dist, norm_dist_transformed, "XYZ", "xyz_hist", "x")
#' }
compare_hist <- function(original, transformed, feature, prefix, xlab) {
  `..count..` <- NULL # Local binding
  ALPHA <- 1
  BINS <- 20
  histogram <- data.frame(
    original = original, 
    transformed = transformed
  )
  original.plot <- ggplot2::ggplot(data = histogram, ggplot2::aes(original)) +
    ggplot2::geom_histogram(alpha = ALPHA, 
                            ggplot2::aes(y = ..count..), 
                            position = "identity", 
                            bins = BINS, 
                            col = "black", 
                            fill = "#FFDF01") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1)) + 
    ggplot2::labs(title = paste("Trait", feature), x = "", y = "") +
    ggplot2::xlab(paste0(feature))
  
  transformed.plot <- ggplot2::ggplot(data = histogram, 
                                      ggplot2::aes(transformed)) +
    ggplot2::geom_histogram(alpha = ALPHA, 
                            ggplot2::aes(y = ..count..), 
                            position = "identity", 
                            bins = BINS, 
                            col = "black", 
                            fill = "#0057A7") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1)) + 
    ggplot2::labs(x="", y="") +
    ggplot2::xlab(latex2exp::TeX(xlab))
  
  ggplot_save(grid::grid.draw(rbind(ggplot2::ggplotGrob(original.plot),
                                    ggplot2::ggplotGrob(transformed.plot), 
                                    size = "last")
                              ),
              paste0(prefix,"_",feature))
}

#' This function generates a single histogram
#'
#' @param data histogram data
#' @param title plot title
#' @param prefix file prefix [default: \code{metapipe}]
#' @param xlab x-axis label [default: \code{NULL}]
#' @param save boolean flag to indicate if that plot should be save to disk [default: \code{TRUE}]
#' @param alpha value to set the transparency of the bins [default: \code{1}]
#' @param angle rotation angle for labels [default: \code{60}]
#' @param bins number of bins to plot [default: \code{20}]
#' @param col colour for bins border [default: \code{"black"}]
#' @param hjust horizontal adjustment of labels [default: \code{1}]
#' @param fill filling colour [default: \code{"#7FCDBB"}]
#' @param is_trait if TRUE then prepends "Trait" to title [default: \code{FALSE}]
#'
#' @export
#'
#' @examples
#' \dontrun{
#'     norm_dist <- rnorm(100)
#'     generate_hist(norm_dist, "XYZ", "xyz_hist", "x")
#'     generate_hist(norm_dist, "XYZ", "xyz_hist", "x", is_trait = TRUE)
#' }
generate_hist <- function(data, 
                          title, 
                          prefix = "metapipe", 
                          xlab = NULL, 
                          save = TRUE, 
                          alpha = 1,
                          angle = 60,
                          bins = 20, 
                          col = "black",
                          hjust = 1,
                          fill = "#7FCDBB",
                          is_trait = FALSE) {
  original <- `..count..` <- NULL
  histogram <- data.frame(original = data)
  # Create tmp title if title != NULL and is_trait = TRUE
  tmp_title <- title
  if (!is.null(title)) {
    tmp_title <- ifelse(is_trait, paste("Trait", title), title)
  }
  
  myPlot <- ggplot2::ggplot(data = histogram, ggplot2::aes(original)) +
    ggplot2::geom_histogram(alpha = alpha, 
                            ggplot2::aes(y = ..count..), 
                            position = "identity", 
                            bins = bins, 
                            col = col, 
                            fill = fill) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angle, 
                                                       hjust = hjust)) + 
    ggplot2::labs(title = tmp_title, x = xlab, y = NULL)

  if (!save)
    return(myPlot)
  filename <- gsub("[[:punct:]]", "", title)
  filename <-  paste0(prefix, "_", gsub(" ", "-", filename))
  ggplot_save(myPlot, filename)
  return(NULL)
}

#' Create Hexagonal logo for MetaPipe
#'
#' @param subplot subplot/image/icon to be used in the background
#' @param dpi dots per inch
#' @param h_color hexagon edge colour
#' @param h_fill hexagon filling colour
#' @param output output filename and path
#' @param p_color package name colour
#'
#' @export
#'
#' @examples
#' \dontrun{
#'     hex_logo()
#' }
hex_logo <- function(subplot = system.file("images/lab-2.png", package = "MetaPipe"),
                     dpi = 800,
                     h_color = "#000000",
                     h_fill = "#363b74",
                     output = system.file("images/metapipe.png", package = "MetaPipe"),
                     p_color = "#eeeeee") {
  hexSticker::sticker(subplot = subplot, package = "MetaPipe", 
                      h_color = h_color,  h_fill = h_fill,
                      dpi = dpi,
                      # l_x = 1.0, l_y = 1.0, spotlight = FALSE, 
                      s_x = 1.0, s_y = .85, s_width = .5,
                      p_x = 1.0, p_y = 1.52, p_size = 6, p_color = p_color,
                      url = "https://github.com/villegar/MetaPipe", 
                      u_angle = 30, u_color = p_color, u_size = 1.35,
                      filename = output)
}
