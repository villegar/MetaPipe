meansp <- read.csv("meansp.csv",header = TRUE)

mean.color <- describeBy(meansp, meansp$FruitColor)   #mean by "Color"
white.mean <- as.data.frame(mean.1$White)         #save these files in wd
write.csv(white.mean, "white.mean.csv")
black.mean <- as.data.frame(mean.1$Black)
write.csv(black.mean, "black.mean.csv")

# color analysis 


# LDA with whole data set

#msp=read.csv("pareto.logdata.colorlab.csv", header = TRUE)
msp <- meansp
library(MASS)
pareto.logdata_color2=msp[ , c(3:1319)]
fit <- lda(FruitColor ~ . , data=pareto.logdata_color2)

library(ggplot2)
lda.data <- cbind(pareto.logdata_color2, predict(fit)$x)
ggplot(lda.data, aes(LD1, LD2)) +
  geom_point(aes(color = FruitColor))


lda.data <- cbind(pareto.logdata_color2, predict(fit)$x)
ggplot(lda.data, aes(LD1, LD2)) +
  geom_point(aes(color = FruitColor))+
  stat_ellipse(aes(x=LD1, y=LD2, fill = FruitColor), alpha = 0.2, geom = "polygon")


# LDA with top 200 (100 from white and 100 from black) 

msp1=read.csv("msp_top200coldif.csv", header = TRUE)
library(MASS)
pareto.logdata_color2=msp1[ , c(4:204)]
fit <- lda(FruitColor ~ . , data=pareto.logdata_color2)



lda.data <- cbind(pareto.logdata_color2, predict(fit)$x)
ggplot(lda.data, aes(LD1, LD2)) +
  geom_point(aes(color = FruitColor))


lda.data <- cbind(pareto.logdata_color2, predict(fit)$x)
ggplot(lda.data, aes(LD1, LD2)) +
  geom_point(aes(color = FruitColor))+
  stat_ellipse(aes(x=LD1, y=LD2, fill = FruitColor), alpha = 0.2, geom = "polygon")