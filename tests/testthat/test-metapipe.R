test_that("load raw data works", {
  set.seed(123)
  # Data frame with the raw data: 1 property and 2 features
  example_data <- data.frame(ID = c(1,2,3,4,5), 
                             P1 = c("one", "two", "three", "four", "five"), 
                             F1 = rnorm(5), 
                             F2 = rnorm(5))
  
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
  example_data <- data.frame(ID = c(1,2,3,4,5), 
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