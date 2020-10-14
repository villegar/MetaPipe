test_that("load raw data works", {
  set.seed(123)
  # Data frame with the raw data: 1 property and 2 traits
  example_data <- data.frame(ID = c(1, 2, 3, 4, 5), 
                             P1 = c("one", "two", "three", "four", "five"), 
                             T1 = rnorm(5), 
                             T2 = rnorm(5),
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
  # Data frame with the raw data: 1 property and 2 traits
  example_data <- data.frame(ID = c(1, 2, 3, 4, 5), 
                             P1 = c("one", "two", "three", "four", "five"), 
                             T1 = rnorm(5), 
                             T2 = rnorm(5),
                             T3 = NA)
 
  # Inserting missing values manually
  example_data$T1[2:3] <- NA
  example_data$T2[4] <- NA
  
  # Expected output when replacing by half of the minimum value
  example_data_alt <- example_data[,]
  example_data_alt$T1[2:3] <- min(example_data$T1, na.rm = TRUE)/2
  example_data_alt$T2[4] <- min(example_data$T2, na.rm = TRUE)/2
  example_data_alt$T3 <- NULL
  
  # Testing the function with different parameters
  expect_message(results_1 <- replace_missing(example_data, c(2)))
  expect_message(results_2 <- replace_missing(example_data, 
                                              c(1, 2), 
                                              prop_na =  0.25))
  expect_message(results_3 <- replace_missing(example_data, 
                                              c(1, 2), 
                                              replace_na =  TRUE))
  
  # Comparing results
  expect_equivalent(example_data[, -5], results_1)       # Drop T3
  expect_equivalent(example_data[, -c(3, 5)], results_2) # Drop T1 and T3
  expect_equivalent(example_data_alt, results_3)  
  
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
                             T1 = rnorm(5), 
                             T2 = rnorm(5))
  expected_output <- data.frame(index = rep(c(1, 2), each = 5),
                                trait = rep(c("T1", "T2"), each = 5),
                                values = c(example_data$T1, example_data$T2),
                                flag = "Normal",
                                transf = "",
                                transf_val = NA,
                                stringsAsFactors = FALSE)
  
  # Transforming data for log_2 normalisation
  example_data_exp2 <- example_data[,]
  example_data_exp2$T1 <- 2 ^ example_data$T1
  
  normalised_data <- assess_normality(example_data, c(1, 2))
  normalised_data2 <- assess_normality(example_data_exp2, 
                                       c(1, 2), 
                                       plots_dir = here::here("plots"))
  
  # Check for generated histograms
  filenames <- c(here::here("plots/HIST_1_LOG_2_T1.png"), 
                 here::here("plots/HIST_2_NORM_T2.png"), 
                 "HIST_1_NORM_T1.png", 
                 "HIST_2_NORM_T2.png")
  for (f in filenames) {
    expect_true(file.exists(f))
    expect_false(dir.exists(f))
    expect_gt(file.size(f), 0)
    file.remove(f)
    expect_false(file.exists(f))
  }
  
  # Delete plots directory
  expect_true(dir.exists(here::here("plots")))
  unlink(here::here("plots"), recursive = TRUE)
  expect_false(dir.exists(here::here("plots")))
  
  # Testing for all data sets
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

test_that("qtl mapping scanone works", {
  # Create toy dataset
  excluded_columns <- c(2)
  population <- 5
  seed <- 123
  set.seed(seed)
  setwd(here::here())
  example_data <- data.frame(ID = 1:population,
                             P1 = c("one", "two", "three", "four", "five"),
                             T1 = rnorm(population),
                             T2 = rnorm(population))
  example_data_normalised <- data.frame(index = rep(c(1, 2), each = 5),
                                        trait = rep(c("T1", "T2"), each = 5),
                                        values = c(example_data$T1, 
                                                   example_data$T2),
                                        flag = "Normal",
                                        transf = "",
                                        transf_val = NA,
                                        stringsAsFactors = FALSE)
  output <- assess_normality(example_data, excluded_columns)
  
  # Create and store random genetic map [for testing only]
  genetic_map <- random_map(population = population, seed = seed)
  write.csv(genetic_map, "metapipe_genetic_map.csv", row.names = FALSE)
  expect_true(file.exists("metapipe_genetic_map.csv"))
  
  x <- qtl::read.cross("csvs", here::here(),
                       genfile = "metapipe_genetic_map.csv",
                       phefile = "metapipe_raw_data_norm.csv")
  traits <- colnames(x$pheno)
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
                             T1 = rnorm(population),
                             T2 = rnorm(population))
  example_data_normalised <- data.frame(index = rep(c(1, 2), each = 5),
                                        trait = rep(c("T1", "T2"), each = 5),
                                        values = c(example_data$T1, 
                                                   example_data$T2),
                                        flag = "Normal",
                                        transf = "",
                                        transf_val = NA,
                                        stringsAsFactors = FALSE)
  
  output <- assess_normality(example_data, excluded_columns)
  
  # Create and store random genetic map [for testing only]
  genetic_map <- random_map(population = population, seed = seed)
  write.csv(genetic_map, "metapipe_genetic_map.csv", row.names = FALSE)
  expect_true(file.exists("metapipe_genetic_map.csv"))
  
  x <- qtl::read.cross("csvs", here::here(),
                       genfile = "metapipe_genetic_map.csv",
                       phefile = "metapipe_raw_data_norm.csv")
  traits <- colnames(x$pheno)
  set.seed(seed)
  x <- qtl::jittermap(x)
  x <- qtl::calc.genoprob(x, step = 1, error.prob = 0.001)
  x_qtl_perm <- qtl_perm_test(x, n_perm = 5, model = "normal", method = "hk")
  expect_equal(c(9, 20), dim(x_qtl_perm))

  filenames <- c("LOD-T1.png",
                 "LOD-T2.png")
  for (f in filenames) {
    expect_true(file.exists(f))
    expect_false(dir.exists(f))
    expect_gt(file.size(f), 0)
    file.remove(f)
    expect_false(file.exists(f))
  }
  
  x_qtl_perm_alt <- qtl_perm_test(x_data = x, 
                                  n_perm = 5, 
                                  raw_data_normalised = example_data_normalised, 
                                  model = "normal", 
                                  method = "hk")
  expect_equal(c(9, 20), dim(x_qtl_perm_alt))
  
  for (f in filenames) {
    expect_true(file.exists(f))
    expect_false(dir.exists(f))
    expect_gt(file.size(f), 0)
    file.remove(f)
    expect_false(file.exists(f))
  }
  
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