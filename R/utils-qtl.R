#' Detect pseudo-markers
#' 
#' Verifies whether or not a string is a pseudo-marker, these contain the 
#' substring \code{'loc'}.
#'
#' @param marker String with the marker.
#'
#' @return \code{TRUE} if the given marker is a pseudo-marker, \code{FALSE} 
#' otherwise.
#'
#' @examples
#' MetaPipe:::is_pseudo_marker('c1.loc1')
#' MetaPipe:::is_pseudo_marker('S1_2345')
#' 
#' @noRd
#' @keywords internal
is_pseudo_marker <- function(marker) {
  return(ifelse(grepl("loc", marker), TRUE, FALSE))
}

#' Transform pseudo-marker
#' 
#' Transform a pseudo-marker to a standard marker by looking for the nearest 
#' markers in a genetic map.
#'
#' @param x_data Cross-data frame containing genetic map data and traits.
#' @param marker Pseudo-marker to be transformed.
#' @param chr Pseudo-marker's chromosome.
#' @param pos Pseudo-marker's position.
#'
#' @return A prime marker and its position.
#'
#' @examples
#' # Toy dataset
#' excluded_columns <- c(1, 2)
#' population <- 5
#' seed <- 1
#' example_data <- data.frame(ID = 1:population,
#'                            P1 = c("one", "two", "three", "four", "five"),
#'                            T1 = rnorm(population),
#'                            T2 = rnorm(population))
#' example_data_normalised <- 
#'   data.frame(index = rep(c(1, 2), each = 5),
#'              trait = rep(c("T1", "T2"), each = 5),
#'              values = c(example_data$T1, example_data$T2),
#'              flag = "Normal",
#'              transf = "",
#'              transf_val = NA,
#'              stringsAsFactors = FALSE)
#' 
#' out_prefix <- here::here("metapipe")
#' example_data_normalised_post <- 
#'   MetaPipe::assess_normality_postprocessing(example_data, 
#'                                             excluded_columns, 
#'                                             example_data_normalised,
#'                                             out_prefix = out_prefix)
#' 
#' # Create and store random genetic map (for testing only)
#' genetic_map <- 
#'   MetaPipe:::random_map(population = population, seed = seed)
#' write.csv(genetic_map, 
#'           here::here("metapipe_genetic_map.csv"), 
#'           row.names = FALSE)
#' # Load cross file with genetic map and raw data for normal traits
#' x <- qtl::read.cross(format = "csvs", 
#'                      dir = here::here(),
#'                      genfile = "metapipe_genetic_map.csv",
#'                      phefile = "metapipe_raw_data_norm.csv")
#' set.seed(seed)
#' x <- qtl::jittermap(x)
#' x <- qtl::calc.genoprob(x, step = 1, error.prob = 0.001)
#' MetaPipe:::transform_pseudo_marker(x, 'loc1', 1, 2.0)
#' 
#' @noRd
#' @keywords internal
transform_pseudo_marker <- function(x_data, marker, chr, pos) {
  markerp <- marker
  posp <- pos
  if (is_pseudo_marker(marker)) {
    minfo <- qtl::find.markerpos(x_data, 
                                 qtl::find.marker(x_data, chr = chr, pos = pos))
    markerp <- rownames(minfo)
    posp <- minfo$pos
  }
  return(c(markerp, as.character(posp)))
}

#' Create effect plots
#' 
#' Create effect plots for significant QTLs found with  
#' \code{\link{qtl_perm_test}}.
#'
#' @param x_data_sim Cross-data frame simulated with \code{qtl::sim.geno}.
#' @param qtl_data Significant QTL data.
#' @param cpus Number of CPUs to be used in the computation.
#' @param plots_dir Output directory for plots.
#'
#' @export
#'
#' @examples
#' # Toy dataset
#' excluded_columns <- c(1, 2)
#' population <- 5
#' seed <- 1
#' example_data <- data.frame(ID = 1:population,
#'                            P1 = c("one", "two", "three", "four", "five"),
#'                            T1 = rnorm(population),
#'                            T2 = rnorm(population))
#' example_data_normalised <- 
#'   data.frame(index = rep(c(1, 2), each = 5),
#'              trait = rep(c("T1", "T2"), each = 5),
#'              values = c(example_data$T1, example_data$T2),
#'              flag = "Normal",
#'              transf = "",
#'              transf_val = NA,
#'              stringsAsFactors = FALSE)
#' 
#' out_prefix <- here::here("metapipe")
#' example_data_normalised_post <- 
#'   MetaPipe::assess_normality_postprocessing(example_data, 
#'                                             excluded_columns, 
#'                                             example_data_normalised,
#'                                             out_prefix = out_prefix)
#' 
#' # Create and store random genetic map (for testing only)
#' genetic_map <- 
#'   MetaPipe:::random_map(population = population, seed = seed)
#' write.csv(genetic_map, 
#'           here::here("metapipe_genetic_map.csv"), 
#'           row.names = FALSE)
#' # Load cross file with genetic map and raw data for normal traits
#' x <- qtl::read.cross(format = "csvs", 
#'                      dir = here::here(),
#'                      genfile = "metapipe_genetic_map.csv",
#'                      phefile = "metapipe_raw_data_norm.csv")
#' set.seed(seed)
#' x <- qtl::jittermap(x)
#' x <- qtl::calc.genoprob(x, step = 1, error.prob = 0.001)
#' x_qtl_perm <- 
#'   MetaPipe::qtl_perm_test(x, n_perm = 5, model = "normal", method = "hk")
#' x_sim <- qtl::sim.geno(x)
#' MetaPipe::effect_plots(x_sim, x_qtl_perm)
#' 
#' @seealso \code{\link{qtl_perm_test}}
effect_plots <- function(x_data_sim, qtl_data, cpus = 1, plots_dir = getwd()) {
  i <- NULL # Local binding
  # Start parallel backend
  cl <- parallel::makeCluster(cpus, setup_strategy = "sequential")
  doParallel::registerDoParallel(cl)
  
  # Load binary operator for backend
  `%dopar%` <- foreach::`%dopar%`
  
  # Extract trait names
  traits <- as.character(qtl_data$trait)
  
  # Extract markers
  markers <- as.character(qtl_data$marker)

  plots <- foreach::foreach(i = 1:nrow(qtl_data)) %dopar% {
    if (qtl_data[i, ]$method == "par-scanone") {
      qtl_data[i, ]$transf <-
        ifelse(is.na(qtl_data[i, ]$transf), "", qtl_data[i, ]$transf)
      if (qtl_data[i, ]$transf == "log") {
        ylab <-
          paste0("$\\log_{", qtl_data[i, ]$transf_val, "}(", traits[i], ")$")
      } else if (qtl_data[i, ]$transf == "root") {
        ylab <-
          paste0("$\\sqrt[", qtl_data[i, ]$transf_val, "]{", traits[i], "}$")
      } else if (qtl_data[i, ]$transf == "power") {
        ylab <- paste0("$(", traits[i], ")^", qtl_data[i, ]$transf_val, "$")
      } else {
        ylab <- traits[i]
      }
      save_plot(
        qtl::effectplot(
          x_data_sim,
          pheno.col = traits[i],
          mname1 = markers[i],
          main = NULL,
          ylab = latex2exp::TeX(ylab)
        ),
        paste0(plots_dir, "/EFF-", traits[i], "-", markers[i])
      )
    } else {
      ylab <- traits[i]
      save_plot(
        qtl::effectplot(
          x_data_sim,
          pheno.col = as.character(traits[i]),
          mname1 = markers[i],
          main = NULL,
          ylab = latex2exp::TeX(ylab)
        ),
        paste0(plots_dir, "/EFF-NP-", traits[i], "-", markers[i])
      )
    }
  }
  parallel::stopCluster(cl) # Stop cluster
}
