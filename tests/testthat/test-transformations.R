test_that("check alpha level works", {
  expect_error(check_alpha(0.5), NA)
  expect_error(check_alpha(-0.5), "alpha must be*")
  expect_error(check_alpha(1.5), "alpha must be*")
  expect_error(check_alpha("1.5"), "alpha must be*")
})
