test_that("is pseudo-marker works", {
  expect_true(is_pseudo_marker('c1.loc1'))
  expect_false(is_pseudo_marker('S1_2345'))
})