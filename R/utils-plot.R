#' Save plot as PNG
#' 
#' Save a graphical object to disk in Portable Network Graphics (PNG) format.
#'
#' @param plt_obj Graphical object.
#' @param name Output name with or without path.
#' @param width Width in inches.
#' @param height Height in inches.
#'
#' @examples
#' MetaPipe:::save_plot(hist(rnorm(100), 
#'                      main = "Histogram of Normal Distribution"), 
#'                      "hist")
#' 
#' @keywords internal
#' @noRd
save_plot <- function(plt_obj, name, width = 6, height = 6) {
  grDevices::png(paste0(name, ".png"), 
                 width = width, 
                 height = height, 
                 units = "in", 
                 res = 300, 
                 type = "cairo")
  print(plt_obj)
  grDevices::dev.off()
}

#' Save plot as TIFF
#' 
#' Save a graphical object to disk in Tagged Image File Format (TIFF) format.
#'
#' @param plt_obj Graphical object.
#' @param name Output name with or without path.
#' @param width Width in inches.
#' @param height Height in inches.
#'
#' @examples
#' MetaPipe:::save_plotTIFF(hist(rnorm(100), 
#'                          main = "Histogram of Normal Distribution"), 
#'                          "hist")
#' 
#' @keywords internal
#' @noRd
save_plotTIFF <- function(plt_obj, name, width = 6, height = 6) {
  grDevices::tiff(paste0(name, ".tiff"), 
                  width = width, 
                  height = height, 
                  units = "in", 
                  res = 300, 
                  type = "cairo")
  print(plt_obj)
  grDevices::dev.off()
}

#' Save plot as PDF
#' 
#' Save a graphical object to disk in Portable Document Format (PDF) format.
#'
#' @param plt_obj Graphical object.
#' @param name Output name with or without path.
#' @param width Width in inches.
#' @param height Height in inches.
#'
#' @examples
#' MetaPipe:::save_plotPDF(hist(rnorm(100), 
#'                         main = "Histogram of Normal Distribution"), 
#'                         "hist")
#' 
#' @keywords internal
#' @noRd
save_plotPDF <- function(plt_obj, name, width = 6, height = 6) {
  grDevices::pdf(paste0(name, ".pdf"),
                 width = width,
                 height = height)
  print(plt_obj)
  grDevices::dev.off()
}

#' Save \code{gpplot2} plot as PNG
#' 
#' Save a graphical object generated with \code{ggplot2} to disk in 
#' Portable Network Graphics (PNG) format.
#'
#' @param plt_obj Graphical object.
#' @param name Output name with or without path.
#' @param width Width in inches.
#' @param height Height in inches.
#'
#' @examples
#' plt_obj <- ggplot2::qplot(rnorm(100))
#' MetaPipe:::ggplot_save(plt_obj, "hist")
#' 
#' @keywords internal
#' @noRd
ggplot_save <- function(plt_obj, name, width = 6, height = 6){
  R.devices::suppressGraphics({
    ggplot2::ggsave(
      paste0(name, ".png"),
      plot   = plt_obj,
      device = "png",
      width  = width,
      height = height,
      dpi    = 300,
      type   = "cairo"
    )
  })
}

#' Compare histograms
#' 
#' Compare two histograms from two different sources. For example original data 
#' and transformed data.
#'
#' @param original Original data.
#' @param transformed Transformed data.
#' @param trait Trait name.
#' @param prefix File prefix.
#' @param xlab X-axis label.
#'
#' @examples
#' norm_dist <- rnorm(100)
#' norm_dist_transformed <- norm_dist^2
#' MetaPipe:::compare_hist(norm_dist, 
#'                         norm_dist_transformed, 
#'                         "XYZ", 
#'                         "xyz_hist", 
#'                         "x")
#'
#' @keywords internal
#' @noRd
compare_hist <- function(original, transformed, trait, prefix, xlab) {
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
    ggplot2::labs(title = paste("Trait", trait), x = "", y = "") +
    ggplot2::xlab(paste0(trait))
  
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
              paste0(prefix,"_",trait))
}

#' Generate histogram
#'
#' @param data Histogram data.
#' @param title Plot title.
#' @param prefix File prefix.
#' @param xlab X-axis label.
#' @param save Boolean flag to indicate if that plot should be save to disk.
#' @param alpha Numeric value to set the transparency of the bins.
#' @param angle Rotation angle for labels.
#' @param bins Number of bins to plot.
#' @param col Colour for bins border.
#' @param hjust Horizontal adjustment of labels.
#' @param fill Filling colour.
#' @param is_trait Boolean flag to prepend "Trait" to \code{title}.
#'
#' @examples
#' norm_dist <- rnorm(100)
#' MetaPipe:::generate_hist(norm_dist, "XYZ", "xyz_hist", "x")
#' MetaPipe:::generate_hist(norm_dist, "XYZ", "xyz_hist", "x", is_trait = TRUE)
#' 
#' @keywords internal
#' @noRd
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
  
  plt_obj <- ggplot2::ggplot(data = histogram, ggplot2::aes(original)) +
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
    return(plt_obj)
  filename <- gsub("[[:punct:]]", "", title)
  filename <-  paste0(prefix, "_", gsub(" ", "-", filename))
  ggplot_save(plt_obj, filename)
  return(NULL)
}

#' Hexagonal logo
#' 
#' Create Hexagonal logo for MetaPipe
#'
#' @param subplot Subplot/image/icon to be used in the background.
#' @param dpi Dots per inch.
#' @param h_color Hexagon edge colour.
#' @param h_fill Hexagon filling colour.
#' @param output Output filename and path.
#' @param p_color Package name colour.
#'
#' @examples
#' MetaPipe:::hex_logo()
#' 
#' @keywords internal
#' @noRd
hex_logo <- function(subplot = system.file("images/lab-2.png", 
                                           package = "MetaPipe"),
                     dpi = 800,
                     h_color = "#000000",
                     h_fill = "#363b74",
                     output = system.file("images/metapipe.png", 
                                          package = "MetaPipe"),
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