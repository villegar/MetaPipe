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