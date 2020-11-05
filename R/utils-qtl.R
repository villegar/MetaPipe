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
#'   MetaPipe:::assess_normality_postprocessing(example_data, 
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
#'   MetaPipe:::assess_normality_postprocessing(example_data, 
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

#' Load data
#' 
#' Load data from a CSV file or verify that input is a data frame.
#'
#' @param input Input data. It can be a string with the filename or a 
#'     data frame.
#' @param wdir Working directory.
#' @param contents Contents of \code{input}.
#' @inheritDotParams utils::read.csv -file
#'
#' @return Data frame.
#' 
#' @keywords internal
load_data <- function(input, wdir = here::here(), contents = "raw", ...) {
    # Verify if input is a string with the filename
    if (class(input) == "character") {
      filename <- file.path(wdir, input)
      # Verify the file exists
      if (!file.exists(filename)) {
        stop("\nThe given path to the file with ", contents, 
             " data was not found: \n",
             filename)
      } else if (!grepl("*.csv$", filename)) { # Verify if input is a CSV file
        stop("\nThe file with ", contents, 
             " data is expected to be in CSV format: \n",
             filename)
      }
      return(read.csv(filename, ...)) # Load data from external file
  } else if(class(input) == "data.frame") { # Verify if input is a data frame
    return(input) # Input was a data frame, so just return original structure
  } else {
    stop("\nValid arguments for ", contents, " data are: \n",
         " - Data frame\n",
         " - Filename (CSV format)")
  }
}

#' Read QTL data
#' 
#' Read QTL data from two sources, one containing genotypes and another one
#' phenotypes.
#'
#' @param geno Data frame or string (filename) to file containing genotype 
#'     data. For example a genetic map like \code{\link{father_riparia}}.
#' @param pheno Data frame or string (filename) to file containing phenotype 
#'     data. For example, output from \code{\link{assess_normality}}.
#' @param wdir Working directory.
#' @param quiet Boolean flag to hide status messages.
#' @param ... Optional parameters for 
#'     \code{\link[qtl:read.cross]{qtl::read.cross(...)}} function.
#'
#' @return Object of \code{cross} class for QTL mapping.
#' @export
#'
#' @examples
#' data(father_riparia)
#' data(ionomics)
#' ionomics_rev <- MetaPipe::replace_missing(ionomics, 
#'                                           excluded_columns = c(1, 2),
#'                                           replace_na =  TRUE)
#' ionomics_normalised <- 
#'   MetaPipe::assess_normality(ionomics_rev,
#'                              excluded_columns = c(1, 2),
#'                              out_prefix = "ionomics",
#'                              transf_vals = c(2, exp(1)),
#'                              show_stats = FALSE)
#' x_data <- MetaPipe::read.cross(father_riparia, 
#'                                ionomics_normalised$norm,
#'                                genotypes = c("nn", "np", "--"))
#' 
#' @family QTL mapping functions
read.cross <- function(geno, pheno, wdir = here::here(), quiet = TRUE, ...) {
  # Local binding
  ID <- NULL
  
  geno <- load_data(geno, wdir, "genotype")
  pheno <- load_data(pheno, wdir, "phenotype")
  
  # First column must be called ID
  colnames(geno)[1] <- colnames(pheno)[1] <- "ID"
  
  # First column must be of the same data type
  geno$ID <- suppressWarnings(as.character(geno$ID))
  pheno$ID <- suppressWarnings(as.character(pheno$ID))
  
  # Check both files have the same IDs
  ids <- dplyr::inner_join(geno, 
                           pheno, 
                           by = "ID")$ID

  if (length(ids) == 0) {
    stop("\nZero matching IDs were found, make sure both genotypes and ",
         "phenotypes use the same ID format.")
  }
  
  # Find matching indices for phenotypes
  idx <- pheno$ID %in% ids
  if (sum(!idx) > 0 && !quiet) {
    message(paste0("\nThe observation",
                   ifelse(sum(!idx) > 1, "s ", " "),
                   "with the following ID",
                   ifelse(sum(!idx) > 1, "s: ", ": "),
                   "\n ", paste0(pheno$ID[!idx], collapse = ", "),
                   "\nwill be dropped, as no matching genotypes were found. "))
  }
  
  # Drop not matching entries (IDs)
  geno <- dplyr::filter(geno, ID %in% ids | is.na(ID) | ID == "")
  geno$ID[is.na(geno$ID)] <- ""
  pheno <- dplyr::filter(pheno, idx)
  
  # Check that working directory exists
  if (!dir.exists(wdir)) {
    stop("\nThe given working directory (wdir) does not exist:\n",
         wdir)
  }
  
  # Create temporal files
  # tmp_dir <- tempdir() 
  write.csv(geno, file.path(wdir, ".geno.csv"), row.names = FALSE)
  write.csv(pheno, file.path(wdir, ".pheno.csv"), row.names = FALSE)
  
  # Load cross file
  x <- qtl::read.cross(format = "csvs", 
                       dir = wdir, 
                       ".geno.csv", 
                       ".pheno.csv",
                       ...)
  
  # Delete temp files
  file.remove(file.path(wdir, ".geno.csv"))
  file.remove(file.path(wdir, ".pheno.csv"))
  return(x)
}
