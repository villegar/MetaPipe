## code to prepare `FatherRipariaMap` dataset
FatherRipariaMap_path <- system.file("extdata", 
                                     "FatherRipariaMap.csv", 
                                     package = "MetaPipe", 
                                     mustWork = TRUE)
FatherRipariaMap <- read.csv(FatherRipariaMap_path)
colnames(FatherRipariaMap)[1] <- "ID"
`%>%` <- dplyr::`%>%`
FatherRipariaMap <- FatherRipariaMap %>%
  dplyr::mutate(ID = gsub("160_271_", "E_", sprintf("%3s", ID)))

# idx <- order(colnames(FatherRipariaMap))
# FatherRipariaMap <- FatherRipariaMap[, idx]
# colnames(FatherRipariaMap)[idx[1:100]]

# FatherRipariaMap$ID <- with(FatherRipariaMap,
#                             gsub("160_271_", "E_", 
#                                  paste0(sprintf("%3s", as.character(ID)))
#                             ))
usethis::use_data(FatherRipariaMap, overwrite = TRUE)
