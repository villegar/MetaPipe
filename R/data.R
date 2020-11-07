#' F1 Seedling Ionomics dataset
#' 
#' A dataset containing information of 21 different isotopes and sample weights  
#' for 217 samples.
#' 
#' @format A data frame with 217 rows and 23 variables:
#' \describe{
#'     \item{ID}{Sample ID.}
#'     \item{SampleWeight}{Sample weight.}
#'     \item{Ca44}{Calcium 44 isotope, in ppm.}
#'     \item{K39}{Potassium 39 isotope, in ppm.}
#'     \item{P31}{Phosphorus 31 isotope, in ppm.}
#'     \item{Li7}{Lithium 7 isotope, in ppm.}
#'     \item{B11}{Boron 11 isotope, in ppm.}
#'     \item{Na23}{Sodium 23, in ppm.}
#'     \item{Mg26}{Magnesium 26, in ppm.}
#'     \item{Al27}{Aluminium 27, in ppm.}
#'     \item{S34}{Sulphur 34, in ppm.}
#'     \item{Fe54}{Iron 54, in ppm.}
#'     \item{Mn55}{Manganese 55, in ppm.}
#'     \item{Co59}{Cobalt 59, in ppm.}
#'     \item{Ni60}{Nickel 60, in ppm.}
#'     \item{Cu63}{Copper 63, in ppm.}
#'     \item{Zn66}{Zinc 66, in ppm.}
#'     \item{As75}{Arsenic 75, in ppm.}
#'     \item{Se78}{Selenium 78, in ppm.}
#'     \item{Rb85}{Rubidium 85, in ppm.}
#'     \item{Sr88}{Strontium 88, in ppm.}
#'     \item{Mo98}{Molybdenum 98, in ppm.}
#'     \item{Cd111}{Cadmium 111, in ppm.}
#' }
#' @author Jason P. Londo \email{jason.londo@usda.gov}
#' @usage data(ionomics)
#' @keywords datasets
# @source \url{https://fennell-lab.netlify.app}
"ionomics"

#' Father Riparia Genetic Map dataset
#' 
#' A dataset containing information of a genetic map for a \emph{vitis Riparia} 
#' subject; more specifically, each column after \code{ID}, contains genotypes 
#' within each chromosome at different positions.
#' 
#' @format A data frame with 296 rows and 1116 variables:
#' \itemize{
#'     \item \strong{ID:} Sample ID.
#'     \item \strong{S\code{<chr>}_\code{<pos>} [1115 variables]:} Chromosome 
#'     \code{<chr>}, position \code{<pos>}.    
#'     \itemize{
#'         \item \strong{S1_61235:}  Chromosome 1, position  61235.
#'         \item \strong{S1_929902:} Chromosome 1, position 929902.
#'         \item \strong{S1_666584:} Chromosome 1, position 666584.
#'    }
#' }
#' @author Gaurab Bhattarai \email{gbhattarai@uga.edu} & 
#' Laszlo G Kovacs \email{laszlokovacs@missouristate.edu}
#' @usage data(father_riparia)
#' @keywords datasets
# @source \url{https://fennell-lab.netlify.app}
#' @references 
#' Bhattarai, Gaurab, "Mapping a New Disease Resistance Locus in an F1 Progeny 
#' Derived from Two Grape Wild Relatives" (2019). MSU Graduate Theses. 3366.
#' \url{https://bearworks.missouristate.edu/theses/3366/}
"father_riparia"