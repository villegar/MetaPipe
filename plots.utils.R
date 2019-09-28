savePlot <- function(myPlot,name, width = 6, height = 6) {
  #pdf(paste0(name,".pdf"), width = 6, height = 6, units = 'in', res = 300)
  png(paste0(name,".png"), width = width, height = height, units = 'in', res = 300, type = "cairo")
  print(myPlot)
  dev.off()
}

savePlotTIFF <- function(myPlot,name, width = 6, height = 6) {
  tiff(paste0(name,".tiff"), width = width, height = height, units = 'in', res = 300, type = "cairo")
  print(myPlot)
  dev.off()
}

savePlotPDF <- function(myPlot,name, width = 6, height = 6) {
  pdf(paste0(name,".pdf"), width = width, height = height)
  print(myPlot)
  dev.off()
}

ggplot.save <- function(myPlot,name){
  suppressGraphics({
    ggsave(
      paste0(name,".png"),
      plot   = myPlot,
      device = 'png',
      width  = 6,
      height = 6,
      dpi    = 300,
      type   = "cairo"
    )
  })
}

# This function generates a two histograms to compare two datasets
## original = original data
## transformed = transformed data
## feature = data feature
## name.prefix = prefix for the file name
## xlab = label for x-axis
generate.overlay.histogram <- function(original,transformed,feature,name.prefix,xlab){
  ALPHA <- 1
  BINS <- 20
  histogram <- data.frame(
    original = original, 
    transformed = transformed
  )
  original.plot <- ggplot(data = histogram, aes(histogram$original)) +
    geom_histogram(alpha = ALPHA, aes(y = ..count..), position = 'identity', bins = BINS, col = "black", fill = "#FFDF01") + #"#B2182B") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
    labs(title=paste("Feature",feature), x='', y='') +
    xlab(paste0(feature))
  
  transformed.plot <- ggplot(data = histogram, aes(histogram$transformed)) +
    geom_histogram(alpha = ALPHA, aes(y = ..count..), position = 'identity', bins = BINS, col = "black", fill = "#0057A7") + #"#2166AC") +
    # geom_density(aes(y=..density..)) + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
    labs(x='', y='') +
    xlab(TeX(xlab))
  
  ggplot.save(grid.draw(rbind(ggplotGrob(original.plot),ggplotGrob(transformed.plot), size = "last")),
              paste0(name.prefix,"_",feature))
}

# This function generates a single histogram
## data = input data
## feature = data feature
## name.prefix = prefix for the file name
## xlab = label for x-axis
generate.histogram <- function(data,feature,name.prefix,xlab){
  ALPHA <- 1
  BINS <- 20
  histogram <- data.frame(original = data)
  myPlot <- ggplot(data = histogram, aes(histogram$original)) +
    geom_histogram(alpha = ALPHA, aes(y = ..count..), position = 'identity', bins = BINS, col = "black", fill = "#7FCDBB") + #"#90C978") +
    # geom_density(aes(y=..density..)) + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
    labs(title=paste("Feature",feature), x='', y='') +
    xlab(paste0(feature))
  ggplot.save(myPlot,paste0(name.prefix,"_",feature))
}
