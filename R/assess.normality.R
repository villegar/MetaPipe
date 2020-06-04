# General libraries
library(ggplot2)
library(grid)
library(gridExtra)
library(latex2exp)
library(R.devices)

# Parallel Libraries
library(foreach)
library(doParallel)

# Source util functions
# source("plots-utils.R")
# source("transformations.R")

normality <- function(){
# Resources
cores <- detectCores()
CPUS <- cores[1] - 1

# Global parameters
plots.directory <- "ionomic"
excluded.columns <- c(1,2)
length.excluded.columns <- length(excluded.columns)
output.files.prefix <- "ionomic"
transformation.values <- c(2,exp(1))#,3,4,5,6,7,8,9,10

# Read and clean dataset
meansp <- read.csv("ionomic.data.csv")
meansp.rows <- nrow(meansp)

print("Starting with Normality Assessment")
dir.create(file.path(getwd(), plots.directory), showWarnings = FALSE)
features <- colnames(meansp)
cl <- makeCluster(CPUS, outfile=paste0('./info_parallel.log'))
registerDoParallel(cl)
transformed.meansp <- foreach(i=(1 + length.excluded.columns):ncol(meansp),
                              .combine =rbind,
                              .packages = c("ggplot2","grid","gridExtra","latex2exp","R.devices")) %dopar% {
                                record <- data.frame(
                                  index = i,
                                  feature = features[i],
                                  values = meansp[,i],
                                  flag = "Non-normal",
                                  transf = "",
                                  transf.value = 0
                                )
                                
                                if(sum(is.finite(meansp[,i]), na.rm = TRUE)>2){
                                  suppressWarnings({
                                    pvalue <- shapiro.test(meansp[,i])[[2]]
                                  })
                                  if(pvalue <= 0.05){ # Data must be transformed
                                    record <- transform_data(pvalue,meansp[,i],features[i],i,length.excluded.columns, plots.directory, transformation.values)
                                    if(length(record)){
                                      record$flag <- "Normal"
                                    }
                                    else {
                                      record <- data.frame(
                                        index = i,
                                        feature = features[i],
                                        values = meansp[,i],
                                        flag = "Non-normal",
                                        transf = "",
                                        transf.value = 0
                                      )
                                    }
                                  }
                                  else{ # Normal data
                                    xlab <- features[i]
                                    transformation <- "NORM"
                                    name.prefix <- paste0(plots.directory,"/HIST_",(i-length.excluded.columns),"_",transformation)
                                    generate.histogram(meansp[,i],features[i],name.prefix,xlab)
                                    record$flag <- "Normal"
                                  }
                                }
                                record
                              }

stopCluster(cl) # Stop cluster
print("Done with Normality Assessment")

normal.transformed.meansp <- transformed.meansp[transformed.meansp$flag == "Normal",]
normal.features <- unique(normal.transformed.meansp$feature)
length.normal.features <- length(normal.features)
non.parametric.meansp <- meansp[,-c(normal.features)]
normal.meansp <- data.frame(matrix(vector(), nrow(normal.transformed.meansp)/length.normal.features, length.normal.features,
                                   dimnames=list(c(), normal.features)),
                            stringsAsFactors = F)
for(i in 1:length.normal.features){
  normal.meansp[i] <- subset(normal.transformed.meansp, feature == normal.features[i])$values
}

# Append excluded columns for transformation
pareto.normal.meansp <- cbind(meansp[,excluded.columns],paretoscale(normal.meansp))
pareto.non.parametric.meansp <- cbind(meansp[,excluded.columns],paretoscale(non.parametric.meansp[,-excluded.columns]))
normal.meansp <- cbind(meansp[,excluded.columns],normal.meansp)

colnames(pareto.normal.meansp)[excluded.columns] <- colnames(meansp)[excluded.columns]
colnames(pareto.non.parametric.meansp)[excluded.columns] <- colnames(meansp)[excluded.columns]
colnames(normal.meansp)[excluded.columns] <- colnames(meansp)[excluded.columns]

write.csv(transformed.meansp, file = paste0(output.files.prefix,".transformed.meansp.csv"))
write.csv(normal.meansp, file = paste0(output.files.prefix,".normal.meansp.csv"))
write.csv(non.parametric.meansp, file = paste0(output.files.prefix,".non.parametric.meansp.csv"))
write.csv(pareto.normal.meansp, file = paste0(output.files.prefix,".pareto.normal.meansp.csv"))
write.csv(pareto.non.parametric.meansp, file = paste0(output.files.prefix,".pareto.non.parametric.meansp.csv"))

# Statistics
normal <- nrow(normal.transformed.meansp[normal.transformed.meansp$transf == "",])/meansp.rows
normal.transformed <- nrow(normal.transformed.meansp)/meansp.rows
total <- nrow(transformed.meansp)/meansp.rows #1316
transformations <- unique(transformed.meansp[c("transf","transf.value")])
transformations <- transformations[-1,]
sorting <- order(transformations$transf, decreasing = T)

transformations <- transformations[sorting,]
cat(paste0("Total features (excluding all NAs features): \t", total))
cat(paste0("\nNormal features (without transformation): \t", normal))
cat(paste0("\nNormal features (transformed): \t\t\t", (normal.transformed - normal)))
cat(paste0("\nTotal Normal features: \t\t\t\t", normal.transformed))
cat(paste0("\nNon-parametric features: \t\t\t", (total - normal.transformed),"\n"))

cat(paste0("\nTransformations summary:"))
cat(paste0("\n\tf(x)\tValue \t# Features"))
for(i in 1:nrow(transformations)){
  cat(paste0("\n\t",transformations$transf[i],"\t",transformations$transf.value[i],"\t"))
  tmp <- subset(normal.transformed.meansp, normal.transformed.meansp$transf == transformations$transf[i])
  tmp <- subset(tmp, transf.value == transformations$transf.value[i])
  cat(nrow(tmp)/meansp.rows)
}
cat("\n\n") # Clean output
}
