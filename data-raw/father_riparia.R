## code to prepare `father_riparia` dataset
father_riparia_path <- system.file("extdata", 
                                     "father_riparia.csv", 
                                     package = "MetaPipe", 
                                     mustWork = TRUE)
father_riparia <- read.csv(father_riparia_path)
colnames(father_riparia)[1] <- "ID"
`%>%` <- dplyr::`%>%`
father_riparia <- father_riparia %>%
  dplyr::mutate(ID = gsub("160_271_", "E_", sprintf("%3s", ID)))

# idx <- order(colnames(father_riparia))
# father_riparia <- father_riparia[, idx]
# colnames(father_riparia)[idx[1:100]]

# father_riparia$ID <- with(father_riparia,
#                             gsub("160_271_", "E_", 
#                                  paste0(sprintf("%3s", as.character(ID)))
#                             ))
usethis::use_data(father_riparia, overwrite = TRUE)
