# Source util functions
source("plots.utils.R")
library(RColorBrewer)

scatterhistY <- function(x, xlab = "", ylab = ""){
  zones = matrix(c(1,2), ncol = 2, byrow = TRUE)
  layout(zones, widths = c(4/5,1/5), heights = c(1))
  xhist <- hist(x, plot = FALSE)
  top <- max(c(xhist$counts))
  par(mar = c(3,3,2,0))
  smoothScatter(x, xaxt='n')
  axis(1, at = seq(1:length(x)), labels = rownames(x), cex = 0.75)
  points(x, col = "red")
  par(mar = c(2,0,1,1))
  barplot(xhist$counts, 
          axes = FALSE, xlim = c(0, top), 
          space = 0, horiz = TRUE)
}

OUT.PREFIX <- "metabolomics"
REPLACE.NA <- FALSE
excluded.columns <- 1:3

# Load raw data
meansp <- read.csv(paste0(OUT.PREFIX,".all.meansp.csv"))
if(!REPLACE.NA){
  NA2halfmin <- function(x) suppressWarnings(replace(x, is.na(x), (min(x, na.rm = TRUE)/2)))
  meansp[,-excluded.columns] <- lapply(meansp[,-excluded.columns], NA2halfmin)
}

transformed.meansp <- meansp # No scaling

transformed.meansp$Group <- as.character(transformed.meansp$Group)
transformed.meansp$Group[is.na(transformed.meansp$Group)] <- "Unknown"
## Calculate mean by color
transformed.meansp.diff.by.color <- 
  data.frame(t(aggregate(transformed.meansp[,-excluded.columns], list(transformed.meansp$Group), mean))[-1,])
transformed.meansp.diff.by.color$X1 <- as.numeric(as.character(transformed.meansp.diff.by.color$X1))
transformed.meansp.diff.by.color$X2 <- as.numeric(as.character(transformed.meansp.diff.by.color$X2))
colnames(transformed.meansp.diff.by.color) <- c("black.mean","white.mean","unknown.mean")
transformed.meansp.diff.by.color$mean.diff <- with(transformed.meansp.diff.by.color, black.mean-white.mean)
transformed.meansp.diff.by.color <- transformed.meansp.diff.by.color[order(transformed.meansp.diff.by.color$black.mean, decreasing = T),]
top.100.black <- rownames(transformed.meansp.diff.by.color)[1:100]
transformed.meansp.diff.by.color <- transformed.meansp.diff.by.color[order(transformed.meansp.diff.by.color$white.mean, decreasing = T),]
top.100.white <- rownames(transformed.meansp.diff.by.color)[1:100]

## Whole dataset and Top 200 features LDA
top.200 <- unique(c(top.100.black,top.100.white))
true.qtl <- read.csv(paste0(OUT.PREFIX,".true.qtl.csv"))
true.qtl.features <- unique(as.character(true.qtl$trait))
pca.top200.true.qtl <- true.qtl.features[true.qtl.features %in% top.200]

meansp.top.features <- meansp[,pca.top200.true.qtl]
rownames(meansp.top.features) <- meansp$ID
meansp.top.features <- meansp.top.features[1:(nrow(meansp)-2),]
meansp.top.features <- meansp.top.features[order(as.numeric(rownames(meansp.top.features))),]
for(i in 1:ncol(meansp.top.features)){
  savePlot(scatterhistY(as.matrix(meansp.top.features[i])),
           paste0("Top200.true.qtl.",colnames(meansp.top.features)[i]))
}

meansp.top.features$ID <- as.numeric(rownames(meansp.top.features))
library(reshape)
library(ggplot2)
library(gridExtra)
df <- melt(meansp.top.features, id=c("ID"))
colnames(df) <- c("ID","Feature","Value")
#df$Value <- with(df,(Value-min(Value))/(max(Value)-min(Value)))
cols <- brewer.pal(n = length(pca.top200.true.qtl), name = "YlGnBu")
for(i in 1:length(pca.top200.true.qtl)){
  savePlot(ggplot(df[df$Feature == pca.top200.true.qtl[i],], aes(Feature, Value, fill = Feature)) + 
             geom_violin(trim = FALSE) + coord_flip() +
             geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, fill = "black" ) + 
             scale_fill_brewer(cols[i]) + theme_bw() + ylab("Ion concentration") + xlab("") + theme(legend.position = "none"),
           paste0("Top200.true.qtl.violin.",pca.top200.true.qtl[i]))
}

for(i in 1:length(pca.top200.true.qtl)){
  savePlot(ggplot(df[df$Feature == pca.top200.true.qtl[i],], aes(Feature, Value, fill = Feature)) + 
             geom_violin(trim = FALSE) + coord_flip() +
             geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, fill = "black" ) + 
             scale_fill_manual(values = c("#CBD3EB")) +
             #scale_fill_brewer(cols[i]) + 
             theme_bw() + ylab("") + xlab("") + theme(legend.position = "none"),
           paste0("Top200.true.qtl.violin.",pca.top200.true.qtl[i]))
}

plots <- vector(mode="list", length = 4)
for(i in 1:length(pca.top200.true.qtl)){
  plots[i] <- ggplot(df[df$Feature == pca.top200.true.qtl[i],], aes(Feature, Value, fill = Feature)) + 
               geom_violin(trim = FALSE) + coord_flip() +
               geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, fill = "black" ) + 
               scale_fill_manual(values = c("#CBD3EB")) +
               #scale_fill_brewer(cols[i]) + 
               theme_bw() + ylab("") + xlab("") + theme(legend.position = "none")
}
DOT.SIZE <- 0.4
plot1 <- ggplot(df[df$Feature == pca.top200.true.qtl[1],], aes(Feature, Value, fill = Feature)) + 
  geom_violin(trim = FALSE) + coord_flip() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = DOT.SIZE, fill = "black" ) + 
  scale_fill_manual(values = c("#CBD3EB")) +
  theme_bw() + ylab("") + xlab("") + theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           axis.text.x = element_blank(),
                                           axis.text.y = element_blank(),
                                           axis.ticks = element_blank()) +
  scale_y_discrete(breaks = c(0,max(df[df$Feature == pca.top200.true.qtl[1],"Value"])),
                   labels = c(0,round(max(df[df$Feature == pca.top200.true.qtl[1],"Value"]))))
plot2 <- ggplot(df[df$Feature == pca.top200.true.qtl[2],], aes(Feature, Value, fill = Feature)) + 
  geom_violin(trim = FALSE) + coord_flip() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = DOT.SIZE, fill = "black" ) + 
  scale_fill_manual(values = c("#CBD3EB")) +
  theme_bw() + ylab("") + xlab("") + theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           axis.text.x = element_blank(),
                                           axis.text.y = element_blank(),
                                           axis.ticks = element_blank()) +
  scale_y_discrete(breaks = c(0,max(df[df$Feature == pca.top200.true.qtl[2],"Value"])), 
                   labels = c(0,round(max(df[df$Feature == pca.top200.true.qtl[2],"Value"]))))
plot3 <- ggplot(df[df$Feature == pca.top200.true.qtl[3],], aes(Feature, Value, fill = Feature)) + 
  geom_violin(trim = FALSE) + coord_flip() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = DOT.SIZE, fill = "black" ) + 
  scale_fill_manual(values = c("#CBD3EB")) +
  theme_bw() + ylab("") + xlab("") + theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           axis.text.x = element_blank(),
                                           axis.text.y = element_blank(),
                                           axis.ticks = element_blank()) +
  scale_y_discrete(breaks = c(0,max(df[df$Feature == pca.top200.true.qtl[3],"Value"])), 
                   labels = c(0,round(max(df[df$Feature == pca.top200.true.qtl[3],"Value"])))) +
  margin(0,0,0,0)
plot4 <- ggplot(df[df$Feature == pca.top200.true.qtl[4],], aes(Feature, Value, fill = Feature)) + 
  geom_violin(trim = FALSE) + coord_flip() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = DOT.SIZE, fill = "black" ) + 
  scale_fill_manual(values = c("#CBD3EB")) +
  theme_bw() + ylab("") + xlab("") + theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           axis.text.x = element_blank(),
                                           axis.text.y = element_blank(),
                                           axis.ticks = element_blank()) +
  scale_y_discrete(breaks = c(0,max(df[df$Feature == pca.top200.true.qtl[4],"Value"])), 
                   labels = c(0,round(max(df[df$Feature == pca.top200.true.qtl[4],"Value"]))))
gA <- ggplotGrob(plot1)
gB <- ggplotGrob(plot2)
gC <- ggplotGrob(plot3)
gD <- ggplotGrob(plot4)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB, gC, gD))
lay <- rbind(c(NA,NA,3,3,3,3,3,NA),
             c(1,2,  3,3,3,3,3,4),
             c(NA,NA,3,3,3,3,3,NA))
lay <- rbind(c(NA,NA,3,3,3,3,3,NA),
             c(1,2,  3,3,3,3,3,4),
             c(NA,NA,3,3,3,3,3,NA),
             c(NA,NA,7,7,7,7,7,NA),
             c(5,6,  7,7,7,7,7,8),
             c(NA,NA,7,7,7,7,7,NA),
             c(NA,NA,11,11,11,11,11,NA),
             c(9,10, 11,11,11,11,11,12),
             c(NA,NA,11,11,11,11,11,NA),
             c(NA,NA,15,15,15,15,15,NA),
             c(13,14,15,15,15,15,15,16),
             c(NA,NA,15,15,15,15,15,NA))
grid.arrange(textGrob(pca.top200.true.qtl[1]), textGrob("0"), plot1, textGrob(max(df[df$Feature == pca.top200.true.qtl[1],"Value"])),
             textGrob(pca.top200.true.qtl[2]), textGrob("0"), plot2, textGrob(max(df[df$Feature == pca.top200.true.qtl[2],"Value"])),
             textGrob(pca.top200.true.qtl[3]), textGrob("0"), plot3, textGrob(max(df[df$Feature == pca.top200.true.qtl[3],"Value"])),
             textGrob(pca.top200.true.qtl[4]), textGrob("0"), plot4, textGrob(max(df[df$Feature == pca.top200.true.qtl[4],"Value"])),
             layout_matrix = lay)

savePlot(grid.arrange(textGrob(pca.top200.true.qtl[1]), textGrob("0"), plot1, textGrob(max(df[df$Feature == pca.top200.true.qtl[1],"Value"])),
                      textGrob(pca.top200.true.qtl[2]), textGrob("0"), plot2, textGrob(max(df[df$Feature == pca.top200.true.qtl[2],"Value"])),
                      textGrob(pca.top200.true.qtl[3]), textGrob("0"), plot3, textGrob(max(df[df$Feature == pca.top200.true.qtl[3],"Value"])),
                      textGrob(pca.top200.true.qtl[4]), textGrob("0"), plot4, textGrob(max(df[df$Feature == pca.top200.true.qtl[4],"Value"])),
                      layout_matrix = lay),
         paste0("grid"),width = 10)
grid.arrange(plot1, plot2, plot3, plot4, nrow = 4)

i <- 1
#df[df$Feature == pca.top200.true.qtl[i],]
ggplot(df[df$Feature != "M338T1357",], aes(Feature, Value, fill = Feature)) + 
  geom_violin(trim = FALSE) + coord_flip() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.4, fill = "black" ) + 
  scale_fill_brewer(palette="YlGnBu") + theme_bw() + ylab("") + xlab("") + theme(legend.position = "bottom", panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                 panel.grid.minor = element_blank())

cat("Feature\t\t Max \t Min")
for(i in 1:length(pca.top200.true.qtl)){
  subsewt <- df[df$Feature == pca.top200.true.qtl[i],]$Value
  cat(pca.top200.true.qtl[i],"\t", max(subsewt), " \t",min(subsewt),"\n")     #0033A0  
}

grid.arrange(p1, p2, nrow = 1)
