args = commandArgs(trailingOnly=TRUE) # Command line arguments
# Install libraries
#install.packages(c("Amelia","ggplot2","grid","gridExtra","latex2exp","psych","R.devices","tidyverse","VIM"), dependencies=T)
#install.packages(c("doParallel","foreach"), dependencies=T)
#install.packages(c("FactoMineR","factoextra"), dependencies=T)
#install.packages(c("qtl"), dependencies=T)
#install.packages("qtl2", repos="http://rqtl.org/qtl2cran")

# General Libraries
library(Amelia, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(grid, quietly = TRUE)
library(gridExtra, quietly = TRUE)
library(latex2exp, quietly = TRUE)
library(plyr, quietly = TRUE)
library(psych, quietly = TRUE)
library(R.devices, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(VIM, quietly = TRUE)

# Parallel Libraries
library(foreach, quietly = TRUE)
library(doParallel, quietly = TRUE)
#library(snow, quietly = TRUE)

# PCA Libraries
library(FactoMineR, quietly = TRUE)
library(factoextra, quietly = TRUE)

# QTL Analysis
library(qtl, quietly = TRUE)

# Source util functions
source("plots.utils.R")
source("qtl.analysis.R")
source("transformations.R")

# Resources
cores <- detectCores()
CPUS <- cores[1] - 1

if(length(args) < 1){
  PERMUTATIONS <- 1000 # Number of permutations for QTL Analysis
  REPLACE.NA <- TRUE
  PARETO.SCALING <- FALSE
  OUT.PREFIX <- "metabolomics"
  PLOTS.DIR <- "metabolomics"
} else if(length(args) < 2){
  PERMUTATIONS <- as.numeric(args[1]) # Number of permutations for QTL Analysis
  REPLACE.NA <- TRUE
  PARETO.SCALING <- FALSE
  OUT.PREFIX <- "metabolomics"
  PLOTS.DIR <- "metabolomics"
} else if(length(args) < 3){
  PERMUTATIONS <- as.numeric(args[1]) # Number of permutations for QTL Analysis
  REPLACE.NA <- as.logical(args[2])
  PARETO.SCALING <- FALSE
  OUT.PREFIX <- "metabolomics"
  PLOTS.DIR <- "metabolomics"
} else if(length(args) < 4){
  PERMUTATIONS <- as.numeric(args[1]) # Number of permutations for QTL Analysis
  REPLACE.NA <- as.logical(args[2])
  PARETO.SCALING <- as.logical(args[3])
  OUT.PREFIX <- "metabolomics"
  PLOTS.DIR <- "metabolomics"
} else if(length(args) < 5){
  PERMUTATIONS <- as.numeric(args[1]) # Number of permutations for QTL Analysis
  REPLACE.NA <- as.logical(args[2])
  PARETO.SCALING <- as.logical(args[3])
  OUT.PREFIX <- args[4]
  PLOTS.DIR <- "metabolomics"
} else {
  PERMUTATIONS <- as.numeric(args[1]) # Number of permutations for QTL Analysis
  REPLACE.NA <- as.logical(args[2])
  PARETO.SCALING <- as.logical(args[3])
  OUT.PREFIX <- args[4]
  PLOTS.DIR <- args[5]
}

cat(paste0("CMD Parameters: (",PERMUTATIONS,",",REPLACE.NA,",",PARETO.SCALING,",",OUT.PREFIX,",",PLOTS.DIR,")"))
# Global parameters
excluded.columns <- c(1,2,3)
length.excluded.columns <- length(excluded.columns)
transformation.values <- c(2,exp(1))#,3,4,5,6,7,8,9,10
input.filename <- "sp.csv"
SEED <- 20190901 # Seed for QTL Analysis
LOD.THRESHOLD <- 3 # LOD threhold for QTL Analysis
NA.COUNT.THRESHOLD <- 0.5 # Allows 50% of NAs per feature


# Environment configuration
dir.create(file.path(getwd(), PLOTS.DIR), showWarnings = FALSE) # Directory for plots

# Load and Cleaning Data
sp <- read.csv(input.filename)
ncols <- ncol(sp)
meansp <- aggregate(sp[,(length.excluded.columns + 1):ncols],by=list(sp$ID),mean, na.action = na.omit)
colnames(meansp)[1] <- "ID"
meansp <- left_join(sp[,c("ID","Group","Generation")],meansp, by="ID")
meansp <- meansp[!duplicated(meansp$ID),]
rownames(meansp) <- 1:nrow(meansp)
meansp.rows <- nrow(meansp)

# Missing Value Plot
#missmap(sp, main = "Missing values vs observed")


# Replacement of Missing Values
# Missing values are replaced by half of the minimum non-zero value for each feature.
if(REPLACE.NA){
  NA2halfmin <- function(x) suppressWarnings(replace(x, is.na(x), (min(x, na.rm = TRUE)/2)))
  meansp[,-excluded.columns] <- lapply(meansp[,-excluded.columns], NA2halfmin)
} else {
  NACount <- which(colMeans(is.na(meansp[,-excluded.columns])) >= NA.COUNT.THRESHOLD) + length.excluded.columns
  if(length(NACount)){
    write.csv(meansp[,c(excluded.columns,NACount)], file = paste0(OUT.PREFIX,".NA.meansp.csv"), row.names=FALSE)
    cat(paste0("The following features were dropped because they have ",(NA.COUNT.THRESHOLD*100),"% or more missing values:\n"))
    cat(colnames(meansp)[NACount])
    meansp[,NACount] <- NULL
  }
}

write.csv(meansp, file = paste0(OUT.PREFIX,".all.meansp.csv"), row.names=FALSE)

# Missing values plot
#missmap(meansp, main = "Missing values vs observed")

generate.boxplots <- function(meansp,ggplot.save){
  print("Generating Boxplots")
  cl <- makeCluster(CPUS, outfile=paste0('./info_parallel.log')) # Make cluster
  registerDoParallel(cl)  # Register cluster
  features <- colnames(meansp)
  AllPlots <- foreach(i=(length.excluded.columns + 1):ncol(meansp), 
                      .packages = c("ggplot2","latex2exp","R.devices")) %dopar% {
                        myPlot <- ggplot(data=meansp,aes(x=ID,y=meansp[,i])) +
                          geom_boxplot(aes(fill= "")) +
                          theme(axis.text.x = element_text(angle = 60, hjust = 1))+ 
                          labs(title=paste("Feature",features[i]), x='ID', y='')
                        ggplot.save(myPlot,paste0("BOX_",(i - length.excluded.columns),"_",features[i]))
                      }
  stopCluster(cl1) # Stop cluster
  print("Done with Boxplots")
}

#generate.boxplots(meansp,ggplot.save)

features <- colnames(meansp)

print("Starting with Normality Assessment")
cl <- makeCluster(CPUS, outfile=paste0('./info_parallel.log'))
registerDoParallel(cl)
transformed.meansp <- foreach(i=(length.excluded.columns + 1):ncol(meansp),
                              .combine =rbind,
                              .packages = c("ggplot2","grid","gridExtra","latex2exp","R.devices")) %dopar% {
                                record <- data.frame( # Create and populate entry for current feature
                                  index = i,
                                  feature = features[i],
                                  values = meansp[,i],
                                  flag = "Non-normal",
                                  transf = "",
                                  transf.value = NA
                                )
                                
                                # Verify the current feature has at least 3 non-NA rows
                                if(sum(is.finite(meansp[,i]), na.rm = TRUE)>2){
                                  pvalue <- shapiro.test(meansp[,i])[[2]] # Assess normality of feature before transforming it
                                  if(pvalue <= 0.05){ # Data must be transformed
                                    record <- transform.data(pvalue,meansp[,i],features[i],i,length.excluded.columns, PLOTS.DIR, transformation.values)
                                    
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
                                        transf.value = NA
                                      )
                                    }
                                  }
                                  else{ # Normal data
                                    xlab <- features[i]
                                    transformation <- "NORM"
                                    name.prefix <- paste0(PLOTS.DIR,"/HIST_",(i - length.excluded.columns),"_",transformation)
                                    generate.histogram(meansp[,i],features[i],name.prefix,xlab)
                                    record$flag <- "Normal"
                                  }
                                }
                                record
                              }

stopCluster(cl) # Stop cluster
print("Done with Normality Assessment")

normal.transformed.meansp <- transformed.meansp[transformed.meansp$flag == "Normal",]
non.parametric.transformed.meansp <- transformed.meansp[transformed.meansp$flag == "Non-normal",]
non.parametric.features <- unique(as.character(non.parametric.transformed.meansp$feature))
normal.features <- unique(as.character(normal.transformed.meansp$feature))
length.normal.features <- length(normal.features)
non.parametric.meansp <- meansp[,non.parametric.features]#meansp[,-c(normal.features)]
normal.meansp <- data.frame(matrix(vector(), nrow(normal.transformed.meansp)/length.normal.features, length.normal.features,
                                   dimnames=list(c(), normal.features)),
                            stringsAsFactors = F)
for(i in 1:length.normal.features){
  normal.meansp[i] <- subset(normal.transformed.meansp, feature == normal.features[i])$values
}

# Append excluded columns for transformation 
if(PARETO.SCALING){ # Apply Pareto Scaling
  transformed.normal.meansp <- cbind(meansp[,excluded.columns],paretoscale(normal.meansp))
  transformed.non.parametric.meansp <- cbind(meansp[,excluded.columns],paretoscale(non.parametric.meansp))
} else { # No Scaling
  transformed.normal.meansp <- cbind(meansp[,excluded.columns],normal.meansp)
  transformed.non.parametric.meansp <- cbind(meansp[,excluded.columns],non.parametric.meansp)
}
normal.meansp <- cbind(meansp[,excluded.columns],normal.meansp)

write.csv(transformed.meansp, file = paste0(OUT.PREFIX,".transformed.all.meansp.csv"), row.names=FALSE)
write.csv(normal.meansp, file = paste0(OUT.PREFIX,".normal.meansp.csv"), row.names=FALSE)
write.csv(non.parametric.meansp, file = paste0(OUT.PREFIX,".non.parametric.meansp.csv"), row.names=FALSE)
write.csv(transformed.normal.meansp, file = paste0(OUT.PREFIX,".transformed.normal.meansp.csv"), row.names=FALSE)
write.csv(transformed.non.parametric.meansp, file = paste0(OUT.PREFIX,".transformed.non.parametric.meansp.csv"), row.names=FALSE)

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

# Prepocessing data for QTL Analysis
geno.map <- read.csv("OriginalMap.csv")
colnames(geno.map)[1] <- "ID"
geno.map$ID <- as.character(geno.map$ID)

## Normal features
transformed.normal.meansp$GenoID <- with(transformed.normal.meansp,
                                         gsub(" ","0",paste0(Generation,"_",sprintf("%3s",as.character(ID))))
)
transformed.normal.meansp$ID <- transformed.normal.meansp$GenoID
transformed.normal.meansp$GenoID <- NULL

normal.phe <- inner_join(transformed.normal.meansp,geno.map, by="ID")[,colnames(transformed.normal.meansp)]
normal.phe$Group <- NULL
normal.phe$Generation <- NULL
normal.gen <- rbind(geno.map[1:2,],inner_join(normal.phe,geno.map, by="ID")[,colnames(geno.map)])


## Non-parametric features
transformed.non.parametric.meansp$GenoID <- with(transformed.non.parametric.meansp,
                                                 gsub(" ","0",paste0(Generation,"_",sprintf("%3s",as.character(ID))))
)
transformed.non.parametric.meansp$ID <- transformed.non.parametric.meansp$GenoID
transformed.non.parametric.meansp$GenoID <- NULL

non.parametric.phe <- inner_join(transformed.non.parametric.meansp,geno.map, by="ID")[,colnames(transformed.non.parametric.meansp)]
non.parametric.phe$Group <- NULL
non.parametric.phe$Generation <- NULL
non.parametric.gen <- rbind(geno.map[1:2,],inner_join(non.parametric.phe,geno.map, by="ID")[,colnames(geno.map)])

# Clean phenotypic data
non.parametric.empty.features <- sapply(non.parametric.phe, function(x) all(is.na(x)) || all(is.infinite(x)))
normal.empty.features <- sapply(normal.phe, function(x) all(is.na(x)) || all(is.infinite(x)))
#non.parametric.phe.ncols <- ncol(non.parametric.phe)
#normal.phe.ncols <- ncol(normal.phe)
non.parametric.phe[non.parametric.empty.features] <- NULL
normal.phe[normal.empty.features] <- NULL

if(any(non.parametric.empty.features)){
  print(paste0("The following non-parametric features were removed (NAs):"))
  print(names(non.parametric.empty.features)[non.parametric.empty.features])
}

if(any(normal.empty.features)){
  print(paste0("The following normal features were removed (NAs):"))
  print(names(normal.empty.features)[normal.empty.features])
}

# Write genotypic and phenotypic dataset
## Normal features
write.csv(normal.gen, file = paste0(OUT.PREFIX,".normal.gen.csv"), row.names=FALSE)
write.csv(normal.phe, file = paste0(OUT.PREFIX,".normal.phe.csv"), row.names=FALSE)
## Non-parametric features
write.csv(non.parametric.gen, file = paste0(OUT.PREFIX,".non.parametric.gen.csv"), row.names=FALSE)
write.csv(non.parametric.phe, file = paste0(OUT.PREFIX,".non.parametric.phe.csv"), row.names=FALSE)


# QTL Analysis
#metabolomics.normal <- read.cross("csvs",".",
#                                paste0(OUT.PREFIX,".normal.gen.csv"),
#                                paste0(OUT.PREFIX,".normal.phe.csv"))
#metabolomics.non.parametric <- read.cross("csvs",".",
#                                paste0(OUT.PREFIX,".non.parametric.gen.csv"),
#                                paste0(OUT.PREFIX,".non.parametric.phe.csv"))

x.normal <- read.cross("csvs",".",
                       paste0(OUT.PREFIX,".normal.gen.csv"),
                       paste0(OUT.PREFIX,".normal.phe.csv"))
features <- colnames(x.normal$pheno)
set.seed(SEED)
x.normal <- jittermap(x.normal)
x.normal <- calc.genoprob(x.normal, step=1, error.prob=0.001)
individuals.phenotyped <- summary(x.normal)[[2]]

x.non.parametric <- read.cross("csvs",".",
                               paste0(OUT.PREFIX,".non.parametric.gen.csv"),
                               paste0(OUT.PREFIX,".non.parametric.phe.csv"))
features.np <- colnames(x.non.parametric$pheno)
x.non.parametric <- jittermap(x.non.parametric)
x.non.parametric <- calc.genoprob(x.non.parametric, step=1, error.prob=0.001)

print("Starting with QTL Analysis")
print("Starting with Normal QTL Analysis")
individuals.phenotyped.np <- summary(x.non.parametric)[[2]]

# Non-MPI cluster for single run Scanone
cl <- makeCluster(CPUS, outfile=paste0('./info_parallel.log'))
registerDoParallel(cl)
x.normal.scanone <- foreach(i=2:ncol(x.normal$pheno),
                            .combine = cbind,
                            .packages = c("ggplot2","grid","gridExtra","latex2exp","qtl","R.devices")) %dopar% {
                              
                              # Run single scan
                              normal.scanone <-  scanone(x.normal, pheno.col = i,  model = "normal", method = "hk")
                              if(i == 2){
                                record <- data.frame(
                                  chr = normal.scanone$chr,
                                  pos = normal.scanone$pos,
                                  lod = normal.scanone$lod,
                                  row.names = rownames(normal.scanone)
                                )
                                colnames(record)[3] <- features[i]
                              }
                              else{
                                record <- data.frame(data = normal.scanone$lod)
                                colnames(record) <- features[i]
                              }
                              record
                            }

x.non.parametric.scanone <- foreach(i=2:ncol(x.non.parametric$pheno),
                                    .combine = cbind,
                                    .packages = c("ggplot2","grid","gridExtra","latex2exp","qtl","R.devices")) %dopar% {
                                      
                                      # Run single scan
                                      non.parametric.scanone <-  scanone(x.non.parametric, pheno.col = i,  model = "np")
                                      if(i == 2){
                                        record <- data.frame(
                                          chr = non.parametric.scanone$chr,
                                          pos = non.parametric.scanone$pos,
                                          lod = non.parametric.scanone$lod,
                                          row.names = rownames(non.parametric.scanone)
                                        )
                                        colnames(record)[3] <- features.np[i]
                                      }
                                      else{
                                        record <- data.frame(data = non.parametric.scanone$lod)
                                        colnames(record) <- features.np[i]
                                      }
                                      record
                                    }
stopCluster(cl) # Stop cluster


library(Rmpi) # for mpi.*
library(snow) # for clusterExport, clusterApply

slaves <- as.numeric(Sys.getenv(c("SLURM_NTASKS")))-1
{ sink("/dev/null"); cl <- makeCluster(slaves, type="MPI"); sink(); } # MPI tasks to use


# Non-parametric QTL
print("Starting with Non-Parametric QTL Analysis")
clusterExport(cl, list('features','individuals.phenotyped','normal.transformed.meansp','x.normal','LOD.THRESHOLD','PERMUTATIONS','PLOTS.DIR','savePlot','x.normal.summary.mapping'))

tick <- proc.time()
result <- clusterApply(cl, 2:max_loop, function(i) x.normal.summary.mapping(i))
x.normal.summary.mapping.values <- do.call(rbind,result)
tock <- proc.time() - tick

stopCluster(cl)
mpi.quit()