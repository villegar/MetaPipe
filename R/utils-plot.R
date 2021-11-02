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
#' MetaPipe:::save_plot(hist(rnorm(100)), 
#'                      file.path(tempdir(), "hist"))
#' 
#' @keywords internal
#' @noRd
save_plot <- function(plt_obj, name, width = 6, height = 6) {
  if (!dir.exists(dirname(name)))
    try(dir.create(dirname(name), recursive = TRUE, showWarnings = FALSE))
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
#' MetaPipe:::save_plotTIFF(hist(rnorm(100)), 
#'                          file.path(tempdir(), "hist"))
#' 
#' @keywords internal
#' @noRd
save_plotTIFF <- function(plt_obj, name, width = 6, height = 6) {
  if (!dir.exists(dirname(name)))
    try(dir.create(dirname(name), recursive = TRUE, showWarnings = FALSE))
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
#' MetaPipe:::save_plotPDF(hist(rnorm(100)), 
#'                         file.path(tempdir(), "hist"))
#' 
#' @keywords internal
#' @noRd
save_plotPDF <- function(plt_obj, name, width = 6, height = 6) {
  if (!dir.exists(dirname(name)))
    try(dir.create(dirname(name), recursive = TRUE, showWarnings = FALSE))
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
#' MetaPipe:::ggplot_save(plt_obj, file.path(tempdir(), "hist"))
#' 
#' @keywords internal
#' @noRd
ggplot_save <- function(plt_obj, name, width = 6, height = 6) {
  if (!dir.exists(dirname(name)))
    try(dir.create(dirname(name), recursive = TRUE, showWarnings = FALSE))
  R.devices::suppressGraphics({
    ggplot2::ggsave(
      filename = paste0(basename(name), ".png"),
      plot   = plt_obj,
      device = "png",
      path   = dirname(name),
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
#' \donttest{
#' norm_dist <- rnorm(100)
#' norm_dist_transformed <- norm_dist^2
#' MetaPipe:::compare_hist(norm_dist, 
#'                         norm_dist_transformed, 
#'                         "XYZ", 
#'                         file.path(tempdir(), "xyz_hist"), 
#'                         "x")
#' }
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
  
  ggplot_save(gridExtra::grid.arrange(original.plot,
                                      transformed.plot,
                                      nrow = 2),
              # grid::grid.draw(rbind(ggplot2::ggplotGrob(original.plot),
              #                       ggplot2::ggplotGrob(transformed.plot), 
              #                       size = "last")
              #                 ),
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
#' \donttest{
#' norm_dist <- rnorm(100)
#' prefix <- file.path(tempdir(), "xyz_hist")
#' MetaPipe:::generate_hist(norm_dist, "XYZ", prefix, "x")
#' MetaPipe:::generate_hist(norm_dist, "XYZ", prefix, "x", is_trait = TRUE)
#' }
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
  # return(NULL)
}

#' Simple PCA
#' 
#' Perform a simple PCA using \code{\link[stats:prcomp]{stats::prcomp}}. 
#' Optionally, it will create a PCA biplot using 
#' \code{\link[factoextra:fviz_pca_biplot]{factoextra::fviz_pca_biplot}} if 
#' \code{plot = TRUE}.
#'
#' @param data A numeric or complex matrix (or data frame) that will be used to
#'     perform the Principal Components Analysis.
#' @param plot Boolean flag to indicate whether or not to create a PCA biplot.
#' @param name Output name with or without path.
#' @param width Width in inches.
#' @param height Height in inches.
#' @param ... Arguments passed on to 
#'     \code{\link[factoextra:fviz_pca_biplot]{factoextra::fviz_pca_biplot}}.
# @inheritDotParams factoextra::fviz_pca_biplot -X -col.var -gradient.cols 
#     -label -addEllipses -repel
#' @return Data frame with PCA result.
#' @export
#'
#' @examples
#' \donttest{
#' out_prefix <- file.path(tempdir(), "metapipe_PCA")
#' # Toy dataset
#' example_data <- data.frame(ID = c(1,2,3,4,5), 
#'                            P1 = c("one", "two", "three", "four", "five"), 
#'                            T1 = rnorm(5), 
#'                            T2 = rnorm(5))
#' example_data_pca <- MetaPipe::PCA(example_data[, -c(1:2)],
#'                                   name = out_prefix)
#' 
#' # F1 Seedling Ionomics dataset
#' data(ionomics) # Includes some missing data
#' ionomics_rev <- MetaPipe::replace_missing(ionomics, 
#'                                           excluded_columns = c(1, 2),
#'                                           replace_na =  TRUE,
#'                                           out_prefix = out_prefix)
#' ionomics_pca <- MetaPipe::PCA(ionomics_rev[, -c(1:2)],
#'                               name = out_prefix)
#' }
PCA <- function(data, name = "PCA", plot = TRUE, width = 6, height = 6, ...) {
  idx <- check_types(data, quiet = FALSE)
  if (length(idx) > 0)
    data <- data[, -idx]
  res.pca <- stats::prcomp(data, scale = TRUE)
  #fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
  #res.pca$eig
  # Biplot with top 10 features 
  # savePlot(fviz_pca_biplot(res.pca, col.var="contrib",
  #                          gradient.cols = c("green","red","blue"),#"#00AFBB" "#E7B800", "#FC4E07"),
  #                          select.var = list(contrib = 2),
  #                          label="var",addEllipses=TRUE, ellipse.level=0.95, repel = TRUE  # Avoid text overlapping
  # ) + xlim(-10, 10) + ylim (-10, 10),
  # paste0(PLOTS.DIR,"/PCA-biplot-top10"),8,8)
  if (plot) {
    grDevices::png(paste0(name, ".png"), 
                   width = width, 
                   height = height, 
                   units = "in", 
                   res = 300, 
                   type = "cairo")
    
    print(
      factoextra::fviz_pca_biplot(
        res.pca,
        col.var = "contrib",
        # gradient.cols = c("green", "yellow", "blue"),
        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
        select.var = list(contrib = 2),
        label = "var",
        addEllipses = TRUE,
        ellipse.level = 0.95,
        repel = TRUE,  # Avoid text overlapping
        ...
      )
    )
    grDevices::dev.off()
  }
  return(res.pca)
}
