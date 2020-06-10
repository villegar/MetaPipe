test_that("plot in PDF format works", {
  save_plotPDF(hist(rnorm(100), main = "Histogram of Normal Distribution"),
              "hist")
  filename <- "hist.pdf"
  expect_true(file.exists(filename))
  expect_false(dir.exists(filename))
  expect_gt(file.size(filename), 0)
  file.remove(filename)
  expect_false(file.exists(filename))
})

test_that("plot in PNG format works", {
  save_plot(hist(rnorm(100), main = "Histogram of Normal Distribution"),
           "hist")
  filename <- "hist.png"
  expect_true(file.exists(filename))
  expect_false(dir.exists(filename))
  expect_gt(file.size(filename), 0)
  file.remove(filename)
  expect_false(file.exists(filename))
})

test_that("plot in TIFF format works", {
  save_plotTIFF(hist(rnorm(100), main = "Histogram of Normal Distribution"),
               "hist")
  filename <- "hist.tiff"
  expect_true(file.exists(filename))
  expect_false(dir.exists(filename))
  expect_gt(file.size(filename), 0)
  file.remove(filename)
  expect_false(file.exists(filename))
})

test_that("compare histograms works", {
  norm_dist <- rnorm(100)
  norm_dist_transformed <- norm_dist^2
  compare_hist(norm_dist, norm_dist_transformed, "XYZ", "hist", "x")
  filename <- "hist_XYZ.png"
  expect_true(file.exists(filename))
  expect_false(dir.exists(filename))
  expect_gt(file.size(filename), 0)
  file.remove(filename)
  expect_false(file.exists(filename))
})

test_that("generate histogram works", {
  norm_dist <- rnorm(100)
  generate_hist(norm_dist, "XYZ", "hist", "x")
  filename <- "hist_XYZ.png"
  expect_true(file.exists(filename))
  expect_false(dir.exists(filename))
  expect_gt(file.size(filename), 0)
  file.remove(filename)
  expect_false(file.exists(filename))
})