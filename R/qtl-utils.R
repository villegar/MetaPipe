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
#' @examples
#' # Create toy dataset
#' excluded_columns <- c(2)
#' population <- 5
#' seed <- 1
#' example_data <- data.frame(ID = 1:population,
#'                            P1 = c("one", "two", "three", "four", "five"),
#'                            F1 = rnorm(population),
#'                            F2 = rnorm(population))
#' example_data_normalised <- data.frame(index = rep(c(1, 2), each = 5),
#'                                       feature = rep(c("F1", "F2"), each = 5),
#'                                       values = c(example_data$F1, example_data$F2),
#'                                       flag = "Normal",
#'                                       transf = "",
#'                                       transf_val = NA,
#'                                       stringsAsFactors = FALSE)
#' assess_normality_postprocessing(example_data, 
#'                                 excluded_columns, 
#'                                 example_data_normalised, 
#'                                 out_prefix = here::here("metapipe"))
#' 
#' # Create and store random genetic map [for testing only]
#' genetic_map <- random_map(population = population, seed = seed)
#' write.csv(genetic_map, here::here("metapipe_genetic_map.csv"), row.names = FALSE)
#' # Load cross file with genetic map and raw data for normal features
#' x <- qtl::read.cross(format = "csvs", 
#'                      dir = here::here(),
#'                      genfile = "metapipe_genetic_map.csv",
#'                      phefile = "metapipe_raw_data_norm.csv")
#' set.seed(seed)
#' x <- qtl::jittermap(x)
#' x <- qtl::calc.genoprob(x, step = 1, error.prob = 0.001)
#' transform_pseudo_marker(x, 'loc1', 1, 2.0)
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

#' Create effect plots for significant QTLS found with  
#' \code{\link{qtl_perm_test}}.
#'
#' @param x_data_sim Cross-data simulated with \code{qtl::sim.geno}
#' @param qtl_data Significant QTL's data
#' @param cpus number of CPUS to be used
#' @param plots_dir output directory for plots
#'
#' @export
#'
#' @examples
#' # Create toy dataset
#' excluded_columns <- c(2)
#' population <- 5
#' seed <- 1
#' example_data <- data.frame(ID = 1:population,
#'                            P1 = c("one", "two", "three", "four", "five"),
#'                            F1 = rnorm(population),
#'                            F2 = rnorm(population))
#' example_data_normalised <- data.frame(index = rep(c(1, 2), each = 5),
#'                                       feature = rep(c("F1", "F2"), each = 5),
#'                                       values = c(example_data$F1, example_data$F2),
#'                                       flag = "Normal",
#'                                       transf = "",
#'                                       transf_val = NA,
#'                                       stringsAsFactors = FALSE)
#' assess_normality_postprocessing(example_data, excluded_columns, 
#'                                 example_data_normalised, 
#'                                 out_prefix = here::here("metapipe"))
#' 
#' # Create and store random genetic map [for testing only]
#' genetic_map <- random_map(population = population, seed = seed)
#' write.csv(genetic_map, here::here("metapipe_genetic_map.csv"), row.names = FALSE)
#' 
#' # Load cross file with genetic map and raw data for normal features
#' x <- qtl::read.cross(format = "csvs", 
#'                      dir = here::here(),
#'                      genfile = "metapipe_genetic_map.csv",
#'                      phefile = "metapipe_raw_data_norm.csv")
#' set.seed(seed)
#' x <- qtl::jittermap(x)
#' x <- qtl::calc.genoprob(x, step = 1, error.prob = 0.001)
#' x_qtl_perm <- MetaPipe::qtl_perm_test(x, n_perm = 5, model = "normal", method = "hk")
#' x_sim <- qtl::sim.geno(x)
#' effect_plots(x_sim, x_qtl_perm)
#' 
#' @seealso \code{\link{qtl_perm_test}}
effect_plots <- function(x_data_sim, qtl_data, cpus = 1, plots_dir = getwd()) {
  # Start parallel backend
  cl <- parallel::makeCluster(cpus, setup_strategy = "sequential")
  doParallel::registerDoParallel(cl)
  
  # Load binary operator for backend
  `%dopar%` <- foreach::`%dopar%`
  
  # Extract feature names
  features <- as.character(qtl_data$trait)
  
  # Extract markers
  markers <- as.character(qtl_data$marker)

  plots <- foreach::foreach(i = 1:nrow(qtl_data)) %dopar% {
    if(qtl_data[i, ]$method == "par-scanone") {
      qtl_data[i, ]$transf <-
        ifelse(is.na(qtl_data[i, ]$transf), "", qtl_data[i, ]$transf)
      if(qtl_data[i, ]$transf == "log") {
        ylab <-
          paste0("$\\log_{", qtl_data[i, ]$transf_val, "}(", features[i], ")$")
      } else if (qtl_data[i, ]$transf == "root") {
        ylab <-
          paste0("$\\sqrt[", qtl_data[i, ]$transf_val, "]{", features[i], "}$")
      } else if (qtl_data[i, ]$transf == "power") {
        ylab <- paste0("$(", features[i], ")^", qtl_data[i, ]$transf_val, "$")
      } else {
        ylab <- features[i]
      }
      MetaPipe::save_plot(
        qtl::effectplot(
          x_data_sim,
          pheno.col = features[i],
          mname1 = markers[i],
          main = NULL,
          ylab = latex2exp::TeX(ylab)
        ),
        paste0(plots_dir, "/EFF-", features[i], "-", markers[i])
      )
    } else {
      ylab <- features[i]
      MetaPipe::save_plot(
        qtl::effectplot(
          x_data_sim,
          pheno.col = as.character(features[i]),
          mname1 = markers[i],
          main = NULL,
          ylab = latex2exp::TeX(ylab)
        ),
        paste0(plots_dir, "/EFF-NP-", features[i], "-", markers[i])
      )
    }
  }
  parallel::stopCluster(cl) # Stop cluster
}
