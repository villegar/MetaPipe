test_dir <- tempdir()
test_that("check alpha level works", {
  expect_error(MetaPipe:::check_alpha(0.5), NA)
  expect_error(MetaPipe:::check_alpha(-0.5), "alpha must be*")
  expect_error(MetaPipe:::check_alpha(1.5), "alpha must be*")
  expect_error(MetaPipe:::check_alpha("1.5"), "alpha must be*")
  expect_error(MetaPipe:::check_alpha("One"), "alpha must be*")
})

test_that("log transformation works", {
  set.seed(123)
  data <- rnorm(100, 5)
  
  plots_dir <- file.path(tempdir(), "plots")
  plots_prefix <- file.path(plots_dir, "HIST")
  expect_message(result_df <- 
                   MetaPipe::log_transformation(2 ^ data, 
                                                "EXP_2",
                                                plots_prefix = plots_prefix))
  expect_equivalent(data, result_df)
  expect_warning(result_df <- log_transformation(data, "EXP_2"))
  # expect_null(result_df)
  filename <- file.path(plots_dir, "HIST_LOG_2_EXP_2.png")
  expect_true(file.exists(filename))
  expect_false(dir.exists(filename))
  expect_gt(file.size(filename), 0)
  file.remove(filename)
  expect_false(file.exists(filename))
})

test_that("power transformation works", {
  set.seed(123)
  data <- rnorm(100, 5)
  
  plots_dir <- file.path(tempdir(), "plots")
  plots_prefix <- file.path(plots_dir, "HIST")
  expect_message(result_df <- 
                   MetaPipe::power_transformation(sqrt(data), 
                                                  "ROOT_2",
                                                  plots_prefix = plots_prefix))
  expect_equivalent(data, result_df)
  expect_warning(result_df <- power_transformation(data, "ROOT_2"))
  # expect_null(result_df)
  filename <- file.path(plots_dir, "HIST_POW_2_ROOT_2.png")
  expect_true(file.exists(filename))
  expect_false(dir.exists(filename))
  expect_gt(file.size(filename), 0)
  file.remove(filename)
  expect_false(file.exists(filename))
})

test_that("root transformation works", {
  set.seed(123)
  data <- rnorm(100, 5)
  
  plots_dir <- file.path(tempdir(), "plots")
  plots_prefix <- file.path(plots_dir, "HIST")
  expect_message(result_df <- 
                   MetaPipe::root_transformation(data ^ 2, 
                                                 "POW_2",
                                                 plots_prefix = plots_prefix))
  expect_equivalent((data ^ 2) ^ (1/exp(1)), result_df)
  expect_warning(result_df <- root_transformation(data, "POW_2"))
  # expect_null(result_df)
  filename <- file.path(plots_dir, "HIST_ROOT_e_POW_2.png")
  expect_true(file.exists(filename))
  expect_false(dir.exists(filename))
  expect_gt(file.size(filename), 0)
  file.remove(filename)
  expect_false(file.exists(filename))
})

test_that("pareto scale works", {
  set.seed(123)
  data <- matrix(rnorm(100, 5), ncol = 2)
  data_new <- MetaPipe:::paretoscale(data)
  expect_equal(mean(data_new), 0)
})

test_that("data transform works", {
  set.seed(123)
  data <- rnorm(100, 5)
  
  # Expected transformation
  ## Exponential data -> Log transformation
  expected_df_log_2 <- data.frame(
    index = "",
    trait  = "EXP_2",
    values = log(2 ^ data, 2),
    flag = "Normal",
    transf = "log",
    transf_val = 2,
    stringsAsFactors = FALSE
  )
  
  ## Squared data -> Power transformation
  expected_df_power_2 <- data.frame(
    index = "",
    trait = "ROOT_2",
    values = (sqrt(data)) ^ 2,
    flag = "Normal",
    transf = "power",
    transf_val = 2,
    stringsAsFactors = FALSE
  )
  
  ## Powered data -> Root transformation
  expected_df_root_2 <- data.frame(
    index = "",
    trait = "POW_2",
    values = (data ^ 2) ^ (1 / exp(1)),
    flag = "Normal",
    transf = "root",
    transf_val = "e",
    stringsAsFactors = FALSE
  )
  
  plots_dir <- file.path(tempdir(), "plots")
  plots_prefix <- file.path(plots_dir, "HIST")
  expect_message(result_df <- 
                   MetaPipe::transform_data(data,
                                            plots_prefix = plots_prefix))
  expect_null(result_df)
  
  result_df_log_2 <- MetaPipe::transform_data(2 ^ data, 
                                              "EXP_2",
                                              plots_prefix = plots_prefix)
  result_df_power_2 <- MetaPipe::transform_data(sqrt(data), 
                                                "ROOT_2",
                                                plots_prefix = plots_prefix)
  result_df_root_2 <- MetaPipe::transform_data(data ^ 2, 
                                               "POW_2",
                                               plots_prefix = plots_prefix)
  
  expect_identical(expected_df_log_2, result_df_log_2)
  expect_identical(expected_df_power_2, result_df_power_2)
  expect_identical(expected_df_root_2, result_df_root_2)
  
  filenames <- file.path(plots_dir,
                         c("HIST_LOG_2_EXP_2.png", 
                           "HIST_POW_2_ROOT_2.png", 
                           "HIST_ROOT_e_POW_2.png"))
  for (f in filenames) {
    expect_true(file.exists(f))
    expect_false(dir.exists(f))
    expect_gt(file.size(f), 0)
    file.remove(f)
    expect_false(file.exists(f))
  }
})
