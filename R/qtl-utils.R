#' Verifies whether or not a string is a pseudo-marker, contains the substring
#' 'loc'
#'
#' @param marker string with the marker
#'
#' @return TRUE or FALSE
#' @export
#'
#' @examples
#' is_pseudo_marker('c1.loc1')
#' is_pseudo_marker('S1_2345')
is_pseudo_marker <- function(marker) {
  return(ifelse(grepl("loc", marker), TRUE, FALSE))
}

#' Transforms a pseudo-marker to a standar marker by looking for the nearest one
#'
#' @param x_data cross-file containing genetic map data and features
#' @param marker pseudo marker to be transformed
#' @param chr pseudo-marker's chromosome 
#' @param pos pseudo-marker's position
#'
#' @return a prime marker and position
#' @export
#'
# @examples
transform_pseudo_marker <- function(x_data, marker, chr, pos) {
  markerp <- marker
  posp <- pos
  if(MetaPipe::is_pseudo_marker(marker)) {
    minfo <- qtl::find.markerpos(x_data, 
                                 qtl::find.marker(x_data, chr = chr, pos = pos))
    markerp <- rownames(minfo)
    posp <- minfo$pos
  }
  return(c(markerp, as.character(posp)))
}