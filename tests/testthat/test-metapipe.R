test_that("load raw data works", {
  set.seed(123)
  # Data frame with the raw data: 1 property and 2 features
  example_data <- data.frame(ID = c(1, 2, 3, 4, 5), 
                             P1 = c("one", "two", "three", "four", "five"), 
                             F1 = rnorm(5), 
                             F2 = rnorm(5),
                             stringsAsFactors = FALSE)
  
  # Raw data without duplicates
  filename <- "example_data.csv"
  write.csv(example_data, filename, row.names = FALSE)
  expect_true(file.exists(filename))
  expect_false(dir.exists(filename))
  expect_gt(file.size(filename), 0)
  
  # Raw data with duplicated entries (rows)
  filename_dup <- "example_data_dup.csv"
  write.csv(example_data[c(1:5, 1, 2), ], filename_dup, row.names = FALSE)
  expect_true(file.exists(filename_dup))
  expect_false(dir.exists(filename_dup))
  expect_gt(file.size(filename_dup), 0)
  
  # Testing the function for both files
  expect_equal(example_data, load_raw("example_data.csv", c(1, 2)))
  expect_equal(example_data, load_raw("example_data_dup.csv", c(1, 2)))
  
  # Deleting raw data files
  file.remove(filename)
  file.remove(filename_dup)
  
  # Testing that the files were deleted
  expect_false(file.exists(filename))
  expect_false(file.exists(filename_dup))
})

test_that("replace missing data works", {
  set.seed(123)
  # Data frame with the raw data: 1 property and 2 features
  example_data <- data.frame(ID = c(1, 2, 3, 4, 5), 
                             P1 = c("one", "two", "three", "four", "five"), 
                             F1 = rnorm(5), 
                             F2 = rnorm(5))
 
  # Inserting missing values manually
  example_data$F1[2:3] <- NA
  example_data$F2[4] <- NA
  
  # Expected output when replacing by half of the minimum value
  example_data_alt <- example_data[,]
  example_data_alt$F1[2:3] <- min(example_data$F1, na.rm = TRUE)/2
  example_data_alt$F2[4] <- min(example_data$F2, na.rm = TRUE)/2
  
  # Testing the function with different parameters
  results_1 <- replace_missing(example_data, c(2))
  expect_message(results_2 <- replace_missing(example_data, c(1, 2), prop_na =  0.25))
  results_3 <- replace_missing(example_data, c(1, 2), replace_na =  TRUE)
  
  # Comparing results
  expect_equal(example_data, results_1)
  expect_equal(example_data[, colnames(example_data) != "F1"], results_2)
  expect_equal(example_data_alt, results_3)
  
  # Checking for file generated in the test where prop_na =  0.25
  filename <- "metapipe_NA_raw_data.csv"
  expect_true(file.exists(filename))
  expect_false(dir.exists(filename))
  expect_gt(file.size(filename), 0)
  file.remove(filename)
  expect_false(file.exists(filename))
})

test_that("normality assessment works", {
  set.seed(123)
  example_data <- data.frame(ID = c(1, 2, 3, 4, 5), 
                             P1 = c("one", "two", "three", "four", "five"), 
                             F1 = rnorm(5), 
                             F2 = rnorm(5))
  expected_output <- data.frame(index = rep(c(1, 2), each = 5),
                                feature = rep(c("F1", "F2"), each = 5),
                                values = c(example_data$F1, example_data$F2),
                                flag = "Normal",
                                transf = "",
                                transf_val = NA,
                                stringsAsFactors = FALSE)
  
  # Transforming data for log_2 normalisation
  example_data_exp2 <- example_data[,]
  example_data_exp2$F1 <- 2 ^ example_data$F1

  # Expected output for log_2 normalisation
  expected_output_exp2 <- expected_output[,]
  expected_output_exp2$transf <- rep(c("log", ""), each = 5)
  expected_output_exp2$transf_val <- rep(c(2, NA), each = 5)
  
  # Testing for both data sets
  expect_equal(expected_output, assess_normality(example_data, c(1, 2)))
  expect_equal(expected_output_exp2, assess_normality(example_data_exp2, c(1, 2)))
  
  # Check for generated histograms
  filenames <- c("HIST_1_LOG_2_F1.png", "HIST_1_NORM_F1.png", "HIST_2_NORM_F2.png")
  for (f in filenames) {
    expect_true(file.exists(f))
    expect_false(dir.exists(f))
    expect_gt(file.size(f), 0)
    file.remove(f)
    expect_false(file.exists(f))
  }
})

test_that("normality assessment postprocessing works", {
  set.seed(123)
  example_data <- data.frame(ID = c(1, 2, 3, 4, 5), 
                             P1 = c("one", "two", "three", "four", "five"), 
                             F1 = rnorm(5), 
                             F2 = rnorm(5))
  expected_output <- data.frame(index = rep(c(1, 2), each = 5),
                                feature = rep(c("F1", "F2"), each = 5),
                                values = c(example_data$F1, example_data$F2),
                                flag = "Normal",
                                transf = "",
                                transf_val = NA,
                                stringsAsFactors = FALSE)
  
  # Transforming data for log_2 normalisation
  example_data_exp2 <- example_data[,]
  example_data_exp2$F1 <- 2 ^ example_data$F1
  
  # Expected output for log_2 normalisation
  expected_output_exp2 <- expected_output[,]
  expected_output_exp2$transf <- rep(c("log", ""), each = 5)
  expected_output_exp2$transf_val <- rep(c(2, NA), each = 5)
  
  # Adding noise to feature to make it non-parametric
  example_data_non_par <- example_data[,]
  example_data_non_par$F1 <- c(0, 15000, 0, 17, 0)
  
  # Expected output for log_2 normalisation
  expected_output_non_par <- expected_output[,]
  expected_output_non_par$flag <- rep(c("Non-normal", "Normal"), each = 5)
  
  # Test format of normalised data
  expect_error(assess_normality_postprocessing(example_data, 
                                               c(1, 2), 
                                               expected_output[, -1]), 
               "raw_data_normalised must be the output*")
  
  # Testing for all data sets
  assess_normality_postprocessing(example_data, c(1, 2), expected_output)
  filenames <- c("metapipe_normalisation_stats.csv", 
                 "metapipe_raw_data_non_par.csv", 
                 "metapipe_raw_data_norm.csv", 
                 "metapipe_raw_data_normalised_all.csv")
  for (f in filenames) {
    expect_true(file.exists(f))
    expect_false(dir.exists(f))
    if (f != "metapipe_raw_data_non_par.csv")
      expect_gt(file.size(f), 0)
    file.remove(f)
    expect_false(file.exists(f))
  }
  
  assess_normality_postprocessing(example_data_exp2, c(1, 2), expected_output_exp2)
  filenames <- c("metapipe_normalisation_stats.csv", 
                 "metapipe_raw_data_non_par.csv", 
                 "metapipe_raw_data_norm.csv", 
                 "metapipe_raw_data_normalised_all.csv")
  for (f in filenames) {
    expect_true(file.exists(f))
    expect_false(dir.exists(f))
    if (f != "metapipe_raw_data_non_par.csv")
      expect_gt(file.size(f), 0)
    file.remove(f)
    expect_false(file.exists(f))
  }
  
  assess_normality_postprocessing(example_data_non_par, c(1, 2), expected_output_non_par)
  filenames <- c("metapipe_normalisation_stats.csv", 
                 "metapipe_raw_data_non_par.csv", 
                 "metapipe_raw_data_norm.csv", 
                 "metapipe_raw_data_normalised_all.csv")
  for (f in filenames) {
    expect_true(file.exists(f))
    expect_false(dir.exists(f))
    if (f != "metapipe_raw_data_non_par.csv")
      expect_gt(file.size(f), 0)
    file.remove(f)
    expect_false(file.exists(f))
  }
})

test_that("normality assessment statistics work", {
  set.seed(123)
  example_data <- data.frame(ID = c(1, 2, 3, 4, 5), 
                             P1 = c("one", "two", "three", "four", "five"), 
                             F1 = rnorm(5), 
                             F2 = rnorm(5))
  expected_output <- data.frame(index = rep(c(1, 2), each = 5),
                                feature = rep(c("F1", "F2"), each = 5),
                                values = c(example_data$F1, example_data$F2),
                                flag = "Normal",
                                transf = "",
                                transf_val = NA,
                                stringsAsFactors = FALSE)
  
  # Transforming data for log_2 normalisation
  example_data_exp2 <- example_data[,]
  example_data_exp2$F1 <- 2 ^ example_data$F1
  
  # Expected output for log_2 normalisation
  expected_output_exp2 <- expected_output[,]
  expected_output_exp2$transf <- rep(c("log", ""), each = 5)
  expected_output_exp2$transf_val <- rep(c(2, NA), each = 5)
  
  # Testing for non-existing input file
  expect_error(assess_normality_stats(), "The file *")
  
  # Testing for all data sets
  assess_normality_postprocessing(example_data, c(1, 2), expected_output)
  filenames <- c("metapipe_normalisation_stats.csv", 
                 "metapipe_raw_data_non_par.csv", 
                 "metapipe_raw_data_norm.csv", 
                 "metapipe_raw_data_normalised_all.csv")
  for (f in filenames) {
    if (f == "metapipe_normalisation_stats.csv")
      expect_message(assess_normality_stats())
    expect_true(file.exists(f))
    expect_false(dir.exists(f))
    if (f != "metapipe_raw_data_non_par.csv")
      expect_gt(file.size(f), 0)
    file.remove(f)
    expect_false(file.exists(f))
  }
  
  assess_normality_postprocessing(example_data_exp2, c(1, 2), expected_output_exp2)
  filenames <- c("metapipe_normalisation_stats.csv", 
                 "metapipe_raw_data_non_par.csv", 
                 "metapipe_raw_data_norm.csv", 
                 "metapipe_raw_data_normalised_all.csv")
  for (f in filenames) {
    if (f == "metapipe_normalisation_stats.csv")
      expect_message(assess_normality_stats())
    expect_true(file.exists(f))
    expect_false(dir.exists(f))
    if (f != "metapipe_raw_data_non_par.csv")
      expect_gt(file.size(f), 0)
    file.remove(f)
    expect_false(file.exists(f))
  }
})

test_that("random genotypes works", {
  expect_equal(c("A","B","A"), random_genotypes(size = 3, seed = 1))
})

test_that("random map works", {
  expected_map <- data.frame(ID = c("", "", 1, 2),
                             S1_1 = c(1, 1, "A", "B"),
                             S1_2 = c(1, 2, "H", "H"),
                             S2_1 = c(2, 1, "A", "H"),
                             S2_2 = c(2, 2, "A", "H"))
  expect_equal(expected_map, random_map(lg = 1:2, markers = 2, population = 2, seed = 123))
})

test_that("qtl mapping scanone works", {
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
  x_norm_scone <- MetaPipe::qtl_scone(x, 1, model = "normal", method = "hk")
  expect_equal(c(190, 4), dim(x_norm_scone))
  
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

test_that("qtl mapping permutation test with scanone works", {
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
  # x_qtl_perm <- MetaPipe::qtl_perm_test(x, n_perm = 5, model = "normal", method = "hk")
  # expect_equal(c(9, 20), dim(x_qtl_perm))
  # 
  # filenames <- c("LOD-F1.png", 
  #                "LOD-F2.png")
  # for (f in filenames) {
  #   expect_true(file.exists(f))
  #   expect_false(dir.exists(f))
  #   expect_gt(file.size(f), 0)
  #   file.remove(f)
  #   expect_false(file.exists(f))
  # }
  
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