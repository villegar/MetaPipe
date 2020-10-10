## code to prepare the `ionomics` dataset
ionomics_path <- system.file("extdata", 
                             "ionomics.csv", 
                             package = "MetaPipe", 
                             mustWork = TRUE)
ionomics <- MetaPipe::load_raw(ionomics_path)
ionomics <- ionomics[order(ionomics$ID), ]
rownames(ionomics) <- seq_len(nrow(ionomics))
# This were inspected manually
ionomics[ionomics[, 6] > 0.10, 6] <- NA
ionomics[ionomics[, 8] > 50, 8] <- NA
ionomics[ionomics[, 10] > 40, 10] <- NA
ionomics[ionomics[, 12] > 100, 12] <- NA
ionomics[ionomics[, 14] > 0.4, 14] <- NA
ionomics[ionomics[, 15] > 1.5, 15] <- NA
ionomics[abs(ionomics[19]) > 1.5, 19] <- NA
usethis::use_data(ionomics, overwrite = TRUE)
