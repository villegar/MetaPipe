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
  
  # Expected output for log_2 normalisation
  expected_output_exp2 <- expected_output[,]
  expected_output_exp2$transf <- rep(c("log", ""), each = 5)
  expected_output_exp2$transf_val <- rep(c(2, NA), each = 5)
  
  # Testing for both data sets
  expect_equal(expected_output, 
               assess_normality_core(example_data, c(1, 2)))
  expect_equal(expected_output_exp2, 
               assess_normality_core(example_data_exp2, c(1, 2)))
  
  # Check for generated histograms
  filenames <- c("HIST_1_LOG_2_T1.png", 
                 "HIST_1_NORM_T1.png", 
                 "HIST_2_NORM_T2.png")
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
  
  # Expected output for log_2 normalisation
  expected_output_exp2 <- expected_output[,]
  expected_output_exp2$transf <- rep(c("log", ""), each = 5)
  expected_output_exp2$transf_val <- rep(c(2, NA), each = 5)
  
  # Adding noise to trait to make it non-parametric
  example_data_non_par <- example_data[,]
  example_data_non_par$T1 <- c(0, 15000, 0, 17, 0)
  
  # Expected output for log_2 normalisation
  expected_output_non_par <- expected_output[,]
  expected_output_non_par$flag <- rep(c("Skewed", "Normal"), each = 5)
  
  # Test format of normalised data
  expect_error(assess_normality_postprocessing(example_data, 
                                               c(1, 2), 
                                               expected_output[, -1]), 
               "raw_data must be the output*")
  
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
  
  assess_normality_postprocessing(example_data_exp2, 
                                  c(1, 2), 
                                  expected_output_exp2)
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
  
  assess_normality_postprocessing(example_data_non_par, 
                                  c(1, 2), 
                                  expected_output_non_par)
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
  
  assess_normality_postprocessing(example_data_exp2, 
                                  c(1, 2), 
                                  expected_output_exp2)
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
