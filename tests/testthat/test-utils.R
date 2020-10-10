test_that("check data types works", {
  example_data <- data.frame(ID = c(1,2,3,4,5), 
                             P1 = c("one", "two", "three", "four", "five"), 
                             T1 = rnorm(5))
  expect_equal(check_types(example_data), 2)
  expect_equal(check_types(example_data, numeric = FALSE), c(1, 3))
  expect_equal(check_types(example_data, excluded_columns = 1, numeric = FALSE),
               c(1, 3))
  expect_warning(check_types(example_data, quiet = FALSE))
})
