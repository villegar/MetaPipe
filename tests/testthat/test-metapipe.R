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
  expect_equal(example_data, load_raw("example_data.csv", c(2)))
  expect_equal(example_data, load_raw("example_data_dup.csv", c(2)))
  
  # Deleting raw data files
  file.remove(filename)
  file.remove(filename_dup)
  
  # Testing that the files were deleted
  expect_false(file.exists(filename))
  expect_false(file.exists(filename_dup))
})