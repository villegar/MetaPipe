test_that("check alpha level works", {
  expect_error(check_alpha(0.5), NA)
  expect_error(check_alpha(-0.5), "alpha must be*")
  expect_error(check_alpha(1.5), "alpha must be*")
  expect_error(check_alpha("1.5"), "alpha must be*")
})

test_that("log transformation works", {
  set.seed(123)
  data <- rnorm(100)
  expected_df <- data.frame(
    flag = "Normal",
    transf = "log",
    transf.value = 2,
    stringsAsFactors = FALSE
  )
  result_df <- log_transformation(2^data, "EXP_2")
  expect_identical(expected_df, result_df)
  expect_warning(result_df <- log_transformation(data, "EXP_2"))
  expect_null(result_df)
  filename <- "HIST_LOG_2_EXP_2.png"
  expect_true(file.exists(filename))
  expect_false(dir.exists(filename))
  expect_gt(file.size(filename), 0)
  file.remove(filename)
  expect_false(file.exists(filename))
})

test_that("power transformation works", {
  set.seed(123)
  data <- rnorm(100, 5)
  expected_df <- data.frame(
    flag = "Normal",
    transf = "power",
    transf.value = 2,
    stringsAsFactors = FALSE
  )
  result_df <- power_transformation(sqrt(data), "ROOT_2")
  expect_identical(expected_df, result_df)
  expect_warning(result_df <- power_transformation(data, "ROOT_2"))
  expect_null(result_df)
  filename <- "HIST_POW_2_ROOT_2.png"
  expect_true(file.exists(filename))
  expect_false(dir.exists(filename))
  expect_gt(file.size(filename), 0)
  file.remove(filename)
  expect_false(file.exists(filename))
})

test_that("root transformation works", {
  set.seed(123)
  data <- rnorm(100, 5)
  expected_df <- data.frame(
    flag = "Normal",
    transf = "root",
    transf.value = "e",
    stringsAsFactors = FALSE
  )
  result_df <- root_transformation(data^2, "POW_2")
  expect_identical(expected_df, result_df)
  expect_warning(result_df <- root_transformation(data, "POW_2"))
  expect_null(result_df)
  filename <- "HIST_ROOT_e_POW_2.png"
  expect_true(file.exists(filename))
  expect_false(dir.exists(filename))
  expect_gt(file.size(filename), 0)
  file.remove(filename)
  expect_false(file.exists(filename))
})
