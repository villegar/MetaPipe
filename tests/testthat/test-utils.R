test_that("check data types works", {
  example_data <- data.frame(ID = c(1,2,3,4,5), 
                             P1 = c("one", "two", "three", "four", "five"), 
                             T1 = rnorm(5))
  expect_equal(MetaPipe:::check_types(example_data), 2)
  # expect_equal(MetaPipe:::check_types(example_data, numeric = FALSE), c(1, 3))
  expect_equal(MetaPipe:::check_types(example_data, 
                                      excluded_columns = 1, 
                                      numeric = FALSE),
               c(1, 3))
  # expect_warning(MetaPipe:::check_types(example_data, quiet = FALSE))
})

test_that("random genotypes works", {
  expect_equal(c("A","B","A"), MetaPipe:::random_genotypes(size = 3, seed = 1))
})

test_that("random map works", {
  expected_map <- data.frame(ID = c("", "", 1, 2),
                             S1_1 = c(1, 1, "A", "B"),
                             S1_2 = c(1, 2, "H", "H"),
                             S2_1 = c(2, 1, "A", "H"),
                             S2_2 = c(2, 2, "A", "H"))
  expect_equal(expected_map, MetaPipe:::random_map(lg = 1:2, 
                                                   markers = 2, 
                                                   population = 2, 
                                                   seed = 123))
})
