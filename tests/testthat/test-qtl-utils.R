test_that("is pseudo-marker works", {
  expect_true(is_pseudo_marker('c1.loc1'))
  expect_false(is_pseudo_marker('S1_2345'))
})

test_that("transform pseudo-marker works", {
  # Create toy dataset
  excluded_columns <- c(2)
  population <- 5
  seed <- 123
  set.seed(seed)
  setwd(here::here())
  example_data <- data.frame(ID = 1:population,
                             P1 = c("one", "two", "three", "four", "five"),
                             F1 = rnorm(population),
                             F2 = rnorm(population))
  example_data_normalised <- data.frame(index = rep(c(1, 2), each = 5),
                                        feature = rep(c("F1", "F2"), each = 5),
                                        values = c(example_data$F1, example_data$F2),
                                        flag = "Normal",
                                        transf = "",
                                        transf_val = NA,
                                        stringsAsFactors = FALSE)
  assess_normality_postprocessing(example_data, excluded_columns, example_data_normalised)
  
  # Create and store random genetic map [for testing only]
  genetic_map <- random_map(population = population, seed = seed)
  write.csv(genetic_map, "metapipe_genetic_map.csv", row.names = FALSE)
  expect_true(file.exists("metapipe_genetic_map.csv"))
  
  x <- qtl::read.cross("csvs", here::here(),
                       genfile = "metapipe_genetic_map.csv",
                       phefile = "metapipe_raw_data_norm.csv")
  features <- colnames(x$pheno)
  set.seed(seed)
  x <- qtl::jittermap(x)
  x <- qtl::calc.genoprob(x, step = 1, error.prob = 0.001)
  markerp <- transform_pseudo_marker(x, 'loc1', 1, 2.0)
  expect_equal(c('S1_2', '2.000001'), markerp)
  
  # Delete temporary files
  filenames <- c("metapipe_normalisation_stats.csv", 
                 "metapipe_raw_data_non_par.csv", 
                 "metapipe_raw_data_norm.csv", 
                 "metapipe_raw_data_normalised_all.csv", 
                 "metapipe_genetic_map.csv")
  for (f in filenames) {
    file.remove(f)
    expect_false(file.exists(f))
  }
})