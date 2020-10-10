## code to prepare the `ionomics` dataset
ionomics_path <- system.file("extdata", 
                             "ionomics.csv", 
                             package = "MetaPipe", 
                             mustWork = TRUE)
ionomics <- MetaPipe::load_raw(ionomics_path)
usethis::use_data(ionomics, overwrite = TRUE)
