test_that("is pseudo-marker works", {
  expect_true(is_pseudo_marker('c1.loc1'))
  expect_false(is_pseudo_marker('S1_2345'))
})

test_that("transform pseudo-marker works", {
  # Create toy dataset
  excluded_columns <- c(1, 2)
  population <- 5
  seed <- 123
  set.seed(seed)
  setwd(here::here())
  example_data <- data.frame(ID = 1:population,
                             P1 = c("one", "two", "three", "four", "five"),
                             T1 = rnorm(population),
                             T2 = rnorm(population))
  example_data_normalised <- data.frame(index = rep(c(1, 2), each = 5),
                                        trait = rep(c("T1", "T2"), each = 5),
                                        values = c(example_data$T1, example_data$T2),
                                        flag = "Normal",
                                        transf = "",
                                        transf_val = NA,
                                        stringsAsFactors = FALSE)
  output <- assess_normality(example_data, excluded_columns)
  
  # Create and store random genetic map [for testing only]
  genetic_map <- MetaPipe:::random_map(population = population, seed = seed)
  write.csv(genetic_map, "metapipe_genetic_map.csv", row.names = FALSE)
  expect_true(file.exists("metapipe_genetic_map.csv"))
  
  x <- qtl::read.cross("csvs", here::here(),
                       genfile = "metapipe_genetic_map.csv",
                       phefile = "metapipe_raw_data_norm.csv")
  traits <- colnames(x$pheno)
  set.seed(seed)
  x <- qtl::jittermap(x)
  x <- qtl::calc.genoprob(x, step = 1, error.prob = 0.001)
  markerp <- transform_pseudo_marker(x, 'loc1', 1, 2.0)
  expect_equal(c('S1_2', '2.000001'), markerp)
  
  # Delete temporary files
  filenames <- c("metapipe_normalisation_stats.csv", 
                 "metapipe_raw_data_non_par.csv", 
                 "metapipe_raw_data_norm.csv", 
                 "metapipe_raw_data_normalised_all.csv", 
                 "metapipe_genetic_map.csv")
  for (f in filenames) {
    file.remove(f)
    expect_false(file.exists(f))
  }
})

test_that("effect plots work", {
  # Create toy dataset
  excluded_columns <- c(1, 2)
  population <- 5
  seed <- 123
  set.seed(seed)
  setwd(here::here())
  example_data <- data.frame(ID = 1:population,
                             P1 = c("one", "two", "three", "four", "five"),
                             T1 = rnorm(population),
                             T2 = rnorm(population))
  example_data_normalised <- data.frame(index = rep(c(1, 2), each = 5),
                                        trait = rep(c("T1", "T2"), each = 5),
                                        values = c(example_data$T1, example_data$T2),
                                        flag = "Normal",
                                        transf = "",
                                        transf_val = NA,
                                        stringsAsFactors = FALSE)
  output <- assess_normality(example_data, excluded_columns)
  
  # Create and store random genetic map [for testing only]
  genetic_map <- MetaPipe:::random_map(population = population, seed = seed)
  write.csv(genetic_map, "metapipe_genetic_map.csv", row.names = FALSE)
  expect_true(file.exists("metapipe_genetic_map.csv"))
  
  x <- qtl::read.cross("csvs", here::here(),
                       genfile = "metapipe_genetic_map.csv",
                       phefile = "metapipe_raw_data_norm.csv")
  traits <- colnames(x$pheno)
  set.seed(seed)
  x <- qtl::jittermap(x)
  x <- qtl::calc.genoprob(x, step = 1, error.prob = 0.001)
  x_qtl_perm <- qtl_perm_test(x, n_perm = 5, model = "normal", method = "hk")
  x_sim <- qtl::sim.geno(x)
  
  # Modify QTL data to include transformation data [for testing only]
  x_qtl_perm[1, c("transf", "transf_val")] <- c("log", "2")
  x_qtl_perm[2, c("transf", "transf_val")] <- c("root", "2")
  x_qtl_perm[3, c("transf", "transf_val")] <- c("power", "2")
  
  # Modify QTL data to include skewed traits
  x_qtl_perm[4, c("method")] <- "skw-scanone"
  
  effect_plots(x_sim, x_qtl_perm)
  
  filenames <- c("EFF-T1-S6_4.png",
                 "EFF-T1-S7_1.png",
                 "EFF-T1-S8_3.png",
                 "EFF-NP-T1-S10_4.png",
                 "EFF-T2-S2_8.png",
                 "EFF-T2-S4_9.png",
                 "EFF-T2-S6_3.png",
                 "EFF-T2-S7_9.png",
                 "EFF-T2-S9_5.png")
  for (f in filenames) {
    expect_true(file.exists(f))
    expect_false(dir.exists(f))
    expect_gt(file.size(f), 0)
    file.remove(f)
    expect_false(file.exists(f))
  }
  
  # Delete temporary files
  filenames <- c("LOD-T1.png",
                 "LOD-T2.png",
                 "metapipe_normalisation_stats.csv", 
                 "metapipe_raw_data_non_par.csv", 
                 "metapipe_raw_data_norm.csv", 
                 "metapipe_raw_data_normalised_all.csv", 
                 "metapipe_genetic_map.csv")
  for (f in filenames) {
    file.remove(f)
    expect_false(file.exists(f))
  }
})

test_that("read cross file works", {
  # Toy dataset
  excluded_columns <- c(1, 2)
  population <- 5
  seed <- 123
  set.seed(seed)
  example_data <- data.frame(ID = 1:population,
                             P1 = c("one", "two", "three", "four", "five"),
                             T1 = rnorm(population),
                             T2 = rnorm(population))
  
  output <- MetaPipe::assess_normality(example_data, 
                                       excluded_columns, 
                                       show_stats = FALSE)
  
  # Create and store random genetic map [for testing only]
  genetic_map <- MetaPipe:::random_map(population = population, seed = seed)
  x_data <- MetaPipe::read.cross(genetic_map, output$norm)
  
  # Remove genotype
  expect_message(x_data <- MetaPipe::read.cross(genetic_map[-3,], 
                                                output$norm,
                                                quiet = FALSE))
  
  genetic_map[2, 2:3] <- 2 # Alter markers position (Warning expected)
  expect_warning(x_data <- MetaPipe::read.cross(genetic_map, output$norm))
})