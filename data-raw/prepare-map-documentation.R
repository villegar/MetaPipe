col_names <- colnames(FatherRipariaMap)[-1]
col_names <- c("S1_61235", "S1_929902", "S1_666584", "S1_695458", "S1_666669")
labels <- paste0("#'     \\item{", 
                 col_names, 
                 "}{Chromosome ", 
                 gsub("\\_..*", "", gsub("S", "", col_names)),
                 ", position ",
                 gsub("S*..\\_", "", col_names),
                 "}")
labels <- c("#'     \\item{ID}{Sample ID}", labels)
# gsub("S*..\\_", "", "S19_23343775")
# gsub("S*..\\_", "", "S1_23343775")
# gsub("\\_..*", "", "19_23343775")
# gsub("\\_..*", "", "1_23343775")
fileConn <- file("output.txt")
writeLines(labels, fileConn)
close(fileConn)

cat(paste0(labels[1:500], collapse = "\n"))
cat(paste0(labels[-c(1:500)], collapse = "\n"))
    