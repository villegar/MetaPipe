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



transform_pseudomarker <- function(cross, marker, chr, pos) {
  markerp <- marker
  posp <- pos
  if(MetaPipe::is_pseudo_marker(marker)) {
    minfo <- qtl::find.markerpos(cross, 
                                 qtl::find.marker(cross, chr = chr, pos = pos))
    markerp <- rownames(minfo)
    posp <- minfo$pos
  }
  return(c(markerp, as.character(posp)))
}