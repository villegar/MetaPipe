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

test_that("plot in PNG format generated with GGPLOT2 works", {
  myplot <- ggplot2::qplot(rnorm(100))
  ggplot_save(myplot, "hist")
  filename <- "hist.png"
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
  generate_hist(norm_dist, "ABC", "hist", "x", is_trait = TRUE)
  output <- generate_hist(norm_dist, "XYZ", "hist", "x", FALSE)
  expect_equal(class(output), c("gg", "ggplot"))
  filename <- c("hist_XYZ.png", "hist_ABC.png")
  for(f in filename) {
    expect_true(file.exists(f))
    expect_false(dir.exists(f))
    expect_gt(file.size(f), 0)
    file.remove(f)
    expect_false(file.exists(f))
  }
})

test_that("hex logo works", {
  filename <- "hex-logo.png"
  hex_logo(output = filename)
  expect_true(file.exists(filename))
  expect_false(dir.exists(filename))
  expect_gt(file.size(filename), 0)
  file.remove(filename)
  expect_false(file.exists(filename))
})