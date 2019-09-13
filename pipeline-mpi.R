args = commandArgs(trailingOnly=TRUE) # Command line arguments
# Install libraries
#install.packages(c("Amelia","ggplot2","grid","gridExtra","latex2exp","psych","R.devices","tidyverse","VIM"), dependencies=T)
#install.packages(c("doParallel","foreach"), dependencies=T)
#install.packages(c("FactoMineR","factoextra"), dependencies=T)
#install.packages(c("qtl"), dependencies=T)
#install.packages("qtl2", repos="http://rqtl.org/qtl2cran")

# General Libraries
library(Amelia)
library(ggplot2)
library(grid)
library(gridExtra)
library(latex2exp)
library(plyr)
library(psych)
library(R.devices)
library(tidyverse)
library(VIM)

# Parallel Libraries
library(foreach)
library(doParallel)
#library(snow)

# PCA Libraries
library(FactoMineR)
library(factoextra)

# QTL Analysis
library(qtl)

# Source util functions
source("plots.utils.R")
source("transformations.R")

# Resources
cores <- detectCores()
CPUS <- cores[1] - 1

if(length(args) < 1){
  PERMUTATIONS <- 1000 # Number of permutations for QTL Analysis
  REPLACE.NA <- TRUE
  PARETO.SCALING <- FALSE
} else if(length(args) < 2){
  PERMUTATIONS <- args[1] # Number of permutations for QTL Analysis
  REPLACE.NA <- TRUE
  PARETO.SCALING <- FALSE
} else if(length(args < 3)){
  PERMUTATIONS <- args[1] # Number of permutations for QTL Analysis
  REPLACE.NA <- args[2]
  PARETO.SCALING <- FALSE
} else {
  PERMUTATIONS <- args[1] # Number of permutations for QTL Analysis
  REPLACE.NA <- args[2]
  PARETO.SCALING <- args[3]
}
# Global parameters
plots.directory <- "metabolomics"
excluded.columns <- c(1,2,3)
length.excluded.columns <- length(excluded.columns)
output.files.prefix <- "metabolomics"
transformation.values <- c(2,exp(1))#,3,4,5,6,7,8,9,10
input.filename <- "sp.csv"
SEED <- 20190901 # Seed for QTL Analysis
LOD.THRESHOLD <- 3 # LOD threhold for QTL Analysis


# Environment configuration
dir.create(file.path(getwd(), plots.directory), showWarnings = FALSE) # Directory for plots

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
}

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
                               record <- transform.data(pvalue,meansp[,i],features[i],i,length.excluded.columns, plots.directory, transformation.values)
                               
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
                               name.prefix <- paste0(plots.directory,"/HIST_",(i - length.excluded.columns),"_",transformation)
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

write.csv(transformed.meansp, file = paste0(output.files.prefix,".transformed.all.meansp.csv"), row.names=FALSE)
write.csv(normal.meansp, file = paste0(output.files.prefix,".normal.meansp.csv"), row.names=FALSE)
write.csv(non.parametric.meansp, file = paste0(output.files.prefix,".non.parametric.meansp.csv"), row.names=FALSE)
write.csv(transformed.normal.meansp, file = paste0(output.files.prefix,".transformed.normal.meansp.csv"), row.names=FALSE)
write.csv(transformed.non.parametric.meansp, file = paste0(output.files.prefix,".transformed.non.parametric.meansp.csv"), row.names=FALSE)

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

# PCA analysis with mean (used no missing data) 

#transformed.normal.meansp <- read.csv("transformed.normal.meansp.csv")
#transformed.normal.meansp$X <- NULL
#transformed.normal.meansp <- transformed.normal.meansp[order(as.character(transformed.normal.meansp$ID)),]


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
write.csv(normal.gen, file = paste0(output.files.prefix,".normal.gen.csv"), row.names=FALSE)
write.csv(normal.phe, file = paste0(output.files.prefix,".normal.phe.csv"), row.names=FALSE)
## Non-parametric features
write.csv(non.parametric.gen, file = paste0(output.files.prefix,".non.parametric.gen.csv"), row.names=FALSE)
write.csv(non.parametric.phe, file = paste0(output.files.prefix,".non.parametric.phe.csv"), row.names=FALSE)


# QTL Analysis
#metabolomics.normal <- read.cross("csvs",".",
#                                paste0(output.files.prefix,".normal.gen.csv"),
#                                paste0(output.files.prefix,".normal.phe.csv"))
#metabolomics.non.parametric <- read.cross("csvs",".",
#                                paste0(output.files.prefix,".non.parametric.gen.csv"),
#                                paste0(output.files.prefix,".non.parametric.phe.csv"))

x.normal <- read.cross("csvs",".",
                paste0(output.files.prefix,".normal.gen.csv"),
                paste0(output.files.prefix,".normal.phe.csv"))
features <- colnames(x.normal$pheno)
set.seed(SEED)
x.normal <- jittermap(x.normal)
x.normal <- calc.genoprob(x.normal, step=1, error.prob=0.001)
individuals.phenotyped <- summary(x.normal)[[2]]

x.non.parametric <- read.cross("csvs",".",
                               paste0(output.files.prefix,".non.parametric.gen.csv"),
                               paste0(output.files.prefix,".non.parametric.phe.csv"))
features.np <- colnames(x.non.parametric$pheno)
x.non.parametric <- jittermap(x.non.parametric)
x.non.parametric <- calc.genoprob(x.non.parametric, step=1, error.prob=0.001)

print("Starting with QTL Analysis")
print("Starting with Normal QTL Analysis")
individuals.phenotyped.np <- summary(x.non.parametric)[[2]]

# Load MPI libraries
library(Rmpi)
library(doMPI)

# In case R exits unexpectedly, have it automatically clean up
# resources taken up by Rmpi (slaves, memory, etc...)
.Last <- function(){
   if (is.loaded("mpi_initialize")){
     if (mpi.comm.size(1) > 0){
       print("Please use mpi.close.Rslaves() to close slaves.")
       mpi.close.Rslaves()
     }
     print("Please use mpi.quit() to quit R")
     .Call("mpi_finalize")
   }
}
# Obtain LOD scores for all features and markers
cl <- startMPIcluster()
registerDoMPI(cl)
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

x.normal.summary.mapping <- foreach(i=2:ncol(x.normal$pheno),
                                    .combine = rbind,
                                    .packages = c("ggplot2","grid","gridExtra","latex2exp","qtl","R.devices")) %dopar% {
                                      transformation.info <- normal.transformed.meansp$feature == features[i]
                                      transformation.info <- normal.transformed.meansp[transformation.info,c("transf","transf.value")][1,]
                                      record <- data.frame(
                                        ID = i - 1,
                                        qtl.ID = NA,
                                        trait = features[i],
                                        ind = individuals.phenotyped,
                                        lg = NA,
                                        lod.peak = NA,
                                        pos.peak = NA,
                                        marker = NA,
                                        pos.p95.bay.int = NA,
                                        marker.p95.bay.int = NA,
                                        pvar = NA,
                                        est.add = NA,
                                        est.dom = NA,
                                        p5.lod.thr = NA,
                                        p10.lod.thr = NA,
                                        p.val = NA,
                                        transf = transformation.info$transf,
                                        transf.val = transformation.info$transf.value,
                                        method = "normal-scanone",
                                        p5.qtl = FALSE,
                                        p10.qtl = FALSE
                                      )
                                      
                                      is.pseudo.marker <- function(marker){
                                        if(grepl("loc",marker)){
                                          return(TRUE)
                                        }
                                        return(FALSE)
                                      }
                                      
                                      transform.pseudomarker <- function(cross, marker, chr, pos){
                                        new.marker <- marker
                                        new.pos <- pos
                                        if(is.pseudo.marker(marker)){
                                          marker.info <- find.markerpos(cross, find.marker(cross, chr = chr, pos = pos))
                                          new.marker <- rownames(marker.info)
                                          new.pos <- marker.info$pos
                                        }
                                        return(c(new.marker,as.character(new.pos)))
                                      }
                                      
                                      # Run single scan
                                      normal.scanone <-  scanone(x.normal, pheno.col = i,  model = "normal", method = "hk")
                                      summary.normal.scanone <- summary(normal.scanone, threshold = LOD.THRESHOLD)
                                      lod.count <- nrow(summary.normal.scanone)
                                      if(lod.count){
                                        for(k in 1:lod.count){
                                          if(k > 1){
                                            #new.record <- record[0,] # Create an empty record object
                                            #new.record[1,] <- NA
                                            new.record <- record[1,] # Create copy of record object
                                            #new.record$ID <- NA # Drop the feature ID
                                          }else{
                                            new.record <- record # Copy record structured and data
                                          }
                                          #lod.count <- sum(normal.scanone$lod >= LOD.THRESHOLD)
                                          
                                          #peak.lod <- normal.scanone$lod == max(normal.scanone$lod)
                                          # Extract Peak QTL information
                                          new.record$lg <- summary.normal.scanone[k,"chr"]
                                          new.record$lod.peak <- summary.normal.scanone[k,"lod"]
                                          new.record$pos.peak <- summary.normal.scanone[k,"pos"]
                                          marker <- rownames(summary.normal.scanone)[k]
                                          # Verify if current QTL has a pseudomarker
                                          marker.info <- transform.pseudomarker(x.normal,marker,new.record$lg,new.record$pos.peak)
                                          new.record$marker <- marker.info[1]
                                          new.record$pos.peak <- as.numeric(marker.info[2])
                                          
                                          if(!is.na(new.record$lg)){
                                            new.record$qtl.ID <- with(new.record, sprintf("%s:%s@%f",features[i],lg,pos.peak))
                                          }
                                          
                                          p95.bayesian <- bayesint(normal.scanone, chr = new.record$lg ,expandtomarkers = TRUE, prob = 0.95)
                                          p95.bayesian <- unique(p95.bayesian)
                                          #p95.bayesian <- summary(normal.scanone,  perms=normal.scanone.per, alpha=0.5, pvalues=TRUE)
                                          low.bound <- 1#p95.bayesian$pos == min(p95.bayesian$pos)
                                          upper.bound <- p95.bayesian$pos == max(p95.bayesian$pos)
                                          
                                          p95.bayesian$marker <- NA # Add new column for markers, prevent duplicated row names
                                          # Verify if the Bayesian interval QTLs have pseudomarkers
                                          for(l in 1:nrow(p95.bayesian)){
                                            marker <- rownames(p95.bayesian)[l]
                                            marker.info <- transform.pseudomarker(x.normal,marker,p95.bayesian[l,"chr"],p95.bayesian[l,"pos"])
                                            p95.bayesian[l,"marker"] <- marker.info[1]
                                            p95.bayesian[l,"pos"] <- as.numeric(marker.info[2])
                                          }
                                          new.record$pos.p95.bay.int <- paste0(p95.bayesian[low.bound,"pos"],"-",
                                                                               p95.bayesian[upper.bound,"pos"])
                                          new.record$marker.p95.bay.int <- paste0(p95.bayesian[low.bound,"marker"],"-",
                                                                                  p95.bayesian[upper.bound,"marker"])
                                          #new.record$marker.p95.bay.int <- paste0(rownames(p95.bayesian)[low.bound],"-",
                                          #                                         rownames(p95.bayesian)[upper.bound])
                                          if(k > 1){
                                            record <- rbind(record,new.record)
                                          }else{
                                            record <- new.record
                                          }
                                        }
                                        
                                        #if(lod.count > 0){
                                        #summary(normal.scanone, threshold = 3)
                                        #lod.plot <- plot(normal.scanone, ylab="LOD Score")
                                        #cat(paste0("Scanone: ",i,"\t\tLODs: ",lod.count,"\n"))
                                        normal.scanone.per <- scanone(x.normal, pheno.col = i, model = "normal", method = "hk", n.perm = PERMUTATIONS)
                                        p5 <- summary(normal.scanone.per)[[1]]  #  5% percent
                                        p10 <- summary(normal.scanone.per)[[2]] # 10% percent
                                        
                                        
                                        lod.plot <- savePlot(plot(normal.scanone, ylab="LOD Score") +
                                                               abline(h=p5, lwd=2, lty="solid", col="red") +
                                                               abline(h=p10, lwd=2, lty="solid", col="red"),
                                                             paste0(plots.directory,"/LOD-",features[i]), width = 18)
                                        
                                        record[,]$p5.lod.thr <- p5
                                        record[,]$p10.lod.thr <- p10
                                        
                                        p5.index <- record$lod.peak >= p5
                                        p10.index <- record$lod.peak >= p10
                                        if(!is.na(p5.index) && any(p5.index)){ record[p5.index,]$p5.qtl <- TRUE }
                                        if(!is.na(p10.index)&& any(p10.index)){ record[p10.index,]$p10.qtl <- TRUE }
                                        
                                        
                                        chr <- as.numeric(summary.normal.scanone$chr)
                                        pos <- as.numeric(summary.normal.scanone$pos)
                                        qtl_s <- makeqtl(x.normal, chr, pos, what=c("prob"))
                                        
                                        for(m in 1:length(chr)){
                                          #qtl_s <- makeqtl(x.normal, chr[m], pos[m], what=c("prob"))
                                          #f <- as.formula(paste0("y~",paste0("Q",seq(1:nrow(summary.normal.scanone)), collapse = " + ")))
                                          f <- as.formula(paste0("y~",paste0("Q",m, collapse = " + ")))
                                          fitqtl <- fitqtl(x.normal, pheno.col = i, qtl_s, formula = f , get.ests = TRUE, model = "normal", method="hk")
                                          summary.fitqtl <- summary(fitqtl)
                                          
                                          if(length(summary.fitqtl)){
                                            p.var <- as.numeric(summary.fitqtl[[1]][1,"%var"])
                                            p.value.f <- as.numeric(summary.fitqtl[[1]][,"Pvalue(F)"])[1]
                                            estimates <- as.numeric(summary.fitqtl$ests[,"est"])[-1]
                                            record[m,]$pvar <- p.var
                                            record[m,]$p.val <- p.value.f
                                            record[m,]$est.add <- estimates[1]
                                            record[m,]$est.dom <- estimates[2]
                                          }
                                        }
                                      }
                                      record
                                    }
# Non-parametric QTL
print("Starting with Non-Parametric QTL Analysis")

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


x.non.parametric.summary.mapping <- foreach(i=2:ncol(x.non.parametric$pheno),
                                            .combine = rbind,
                                            .packages = c("ggplot2","grid","gridExtra","latex2exp","qtl","R.devices")) %dopar% {
                                              #transformation.info <- non.parametric.transformed.meansp$feature == features.np[i]
                                              #transformation.info <- non.parametric.transformed.meansp[transformation.info,c("transf","transf.value")][1,]
                                              
                                              record <- data.frame(
                                                ID = i - 1,
                                                qtl.ID = NA,
                                                trait = features.np[i],
                                                ind = individuals.phenotyped.np,
                                                lg = NA,
                                                lod.peak = NA,
                                                pos.peak = NA,
                                                marker = NA,
                                                pos.p95.bay.int = NA,
                                                marker.p95.bay.int = NA,
                                                #pvar = NA,
                                                #est.add = NA,
                                                #est.dom = NA,
                                                p5.lod.thr = NA,
                                                p10.lod.thr = NA,
                                                #p.val = NA,
                                                #transf = transformation.info$transf,
                                                #transf.val = transformation.info$transf.value,
                                                method = "non.parametric-scanone",
                                                p5.qtl = FALSE,
                                                p10.qtl = FALSE
                                              )
                                              
                                              is.pseudo.marker <- function(marker){
                                                if(grepl("loc",marker)){
                                                  return(TRUE)
                                                }
                                                return(FALSE)
                                              }
                                              
                                              transform.pseudomarker <- function(cross, marker, chr, pos){
                                                new.marker <- marker
                                                new.pos <- pos
                                                if(is.pseudo.marker(marker)){
                                                  marker.info <- find.markerpos(cross, find.marker(cross, chr = chr, pos = pos))
                                                  new.marker <- rownames(marker.info)
                                                  new.pos <- marker.info$pos
                                                }
                                                return(c(new.marker,as.character(new.pos)))
                                              }
                                              
                                              # Run single scan
                                              non.parametric.scanone <-  scanone(x.non.parametric, pheno.col = i,  model = "np")
                                              summary.non.parametric.scanone <- summary(non.parametric.scanone, threshold = LOD.THRESHOLD)
                                              lod.count <- nrow(summary.non.parametric.scanone)
                                              if(lod.count){
                                                for(k in 1:lod.count){
                                                  if(k > 1){
                                                    #new.record <- record[0,] # Create an empty record object
                                                    #new.record[1,] <- NA
                                                    new.record <- record[1,] # Create copy of record object
                                                    #new.record$ID <- NA # Drop the feature ID
                                                  }else{
                                                    new.record <- record # Copy record structured and data
                                                  }
                                                  #lod.count <- sum(non.parametric.scanone$lod >= LOD.THRESHOLD)
                                                  
                                                  #peak.lod <- non.parametric.scanone$lod == max(non.parametric.scanone$lod)
                                                  # Extract Peak QTL information
                                                  new.record$lg <- summary.non.parametric.scanone[k,"chr"]
                                                  new.record$lod.peak <- summary.non.parametric.scanone[k,"lod"]
                                                  new.record$pos.peak <- summary.non.parametric.scanone[k,"pos"]
                                                  marker <- rownames(summary.non.parametric.scanone)[k]
                                                  # Verify if current QTL has a pseudomarker
                                                  marker.info <- transform.pseudomarker(x.non.parametric,marker,new.record$lg,new.record$pos.peak)
                                                  new.record$marker <- marker.info[1]
                                                  new.record$pos.peak <- as.numeric(marker.info[2])
                                                  
                                                  if(!is.na(new.record$lg)){
                                                    new.record$qtl.ID <- with(new.record, sprintf("%s:%s@%f",features.np[i],lg,pos.peak))
                                                  }
                                                  
                                                  p95.bayesian <- bayesint(non.parametric.scanone, chr = new.record$lg ,expandtomarkers = TRUE, prob = 0.95)
                                                  p95.bayesian <- unique(p95.bayesian)
                                                  #p95.bayesian <- summary(non.parametric.scanone,  perms=non.parametric.scanone.per, alpha=0.5, pvalues=TRUE)
                                                  low.bound <- 1#p95.bayesian$pos == min(p95.bayesian$pos)
                                                  upper.bound <- p95.bayesian$pos == max(p95.bayesian$pos)
                                                  
                                                  p95.bayesian$marker <- NA # Add new column for markers, prevent duplicated row names
                                                  # Verify if the Bayesian interval QTLs have pseudomarkers
                                                  for(l in 1:nrow(p95.bayesian)){
                                                    marker <- rownames(p95.bayesian)[l]
                                                    marker.info <- transform.pseudomarker(x.non.parametric,marker,p95.bayesian[l,"chr"],p95.bayesian[l,"pos"])
                                                    p95.bayesian[l,"marker"] <- marker.info[1]
                                                    p95.bayesian[l,"pos"] <- as.numeric(marker.info[2])
                                                  }
                                                  new.record$pos.p95.bay.int <- paste0(p95.bayesian[low.bound,"pos"],"-",
                                                                                       p95.bayesian[upper.bound,"pos"])
                                                  new.record$marker.p95.bay.int <- paste0(p95.bayesian[low.bound,"marker"],"-",
                                                                                          p95.bayesian[upper.bound,"marker"])
                                                  
                                                  if(k > 1){
                                                    record <- rbind(record,new.record)
                                                  }else{
                                                    record <- new.record
                                                  }
                                                }
                                                
                                                non.parametric.scanone.per <- scanone(x.non.parametric, pheno.col = i, model = "np", n.perm = PERMUTATIONS)
                                                p5 <- summary(non.parametric.scanone.per)[[1]]  #  5% percent
                                                p10 <- summary(non.parametric.scanone.per)[[2]] # 10% percent
                                                
                                                
                                                lod.plot <- savePlot(plot(non.parametric.scanone, ylab="LOD Score") +
                                                                       abline(h=p5, lwd=2, lty="solid", col="red") +
                                                                       abline(h=p10, lwd=2, lty="solid", col="red"),
                                                                     paste0(plots.directory,"/LOD-NP-",features.np[i]), width = 18)
                                                
                                                record[,]$p5.lod.thr <- p5
                                                record[,]$p10.lod.thr <- p10
                                                
                                                p5.index <- record$lod.peak >= p5
                                                p10.index <- record$lod.peak >= p10
                                                if(!is.na(p5.index) && any(p5.index)){ record[p5.index,]$p5.qtl <- TRUE }
                                                if(!is.na(p10.index)&& any(p10.index)){ record[p10.index,]$p10.qtl <- TRUE }
                                                
                                              }
                                              record
                                            }

# Generate effect plots
x2.non.parametric <- sim.geno(x.non.parametric)
x2.normal <- sim.geno(x.normal)
#effectplot(x2, pheno.col = "M155T28", mname1 = "gbs_13_305342", main = NULL)
t.qtl <- rbind.fill(x.normal.summary.mapping[x.normal.summary.mapping$p5.qtl,],
                    x.non.parametric.summary.mapping[x.non.parametric.summary.mapping$p5.qtl,])
threshold3.qtl <- rbind.fill(x.normal.summary.mapping[x.normal.summary.mapping$lod.peak > LOD.THRESHOLD,],
                             x.non.parametric.summary.mapping[x.non.parametric.summary.mapping$lod.peak > THRESHOLD,])
features.t.qtl <- as.character(t.qtl$trait)
markers.t.qtl <- as.character(t.qtl$marker)

effect.plots <- foreach(i=1:nrow(t.qtl),
        .packages = c("latex2exp","qtl","R.devices")) %dopar% {
          if(t.qtl[i,]$method == "normal-scanone"){
            if(t.qtl[i,]$transf == "log"){
              ylab <- paste0("$\\log_{",t.qtl[i,]$transf.val,"}(",features.t.qtl[i],")$")
            } else if(t.qtl[i,]$transf == "root"){
              ylab <- paste0("$\\sqrt[",t.qtl[i,]$transf.val,"]{",features.t.qtl[i],"}$")
            } else if(t.qtl[i,]$transf == "power"){
              ylab <- paste0("$(",features.t.qtl[i],")^",t.qtl[i,]$transf.val,"$")
            } else {
              ylab <- features.t.qtl[i]
            }
            effect.plot <- savePlot(effectplot(x2.normal, pheno.col = features.t.qtl[i], 
                                               mname1 = markers.t.qtl[i], main = NULL, ylab = TeX(ylab)),
                                    paste0(plots.directory,"/EFF-",features.t.qtl[i],"-",markers.t.qtl[i]))
          } else {
            ylab <- features.t.qtl[i]
            effect.plot <- savePlot(effectplot(x2.non.parametric, pheno.col = as.character(features.t.qtl[i]), 
                                               mname1 = markers.t.qtl[i], main = NULL, ylab = TeX(ylab)),
                                    paste0(plots.directory,"/EFF-NP-",features.t.qtl[i],"-",markers.t.qtl[i]))
          }
        }

closeCluster(cl) # Stop cluster

write.csv(x.normal.scanone, file = paste0(output.files.prefix,".normal.scanone.csv"))
write.csv(x.normal.summary.mapping, file = paste0(output.files.prefix,".normal.summary.mapping.csv"), row.names=FALSE, na="")
write.csv(x.non.parametric.scanone, file = paste0(output.files.prefix,".non.parametric.scanone.csv"))
write.csv(x.non.parametric.summary.mapping, file = paste0(output.files.prefix,".non.parametric.summary.mapping.csv"), row.names=FALSE, na="")
write.csv(t.qtl, file = paste0(output.files.prefix,".true.qtl.csv"), row.names=FALSE, na="")
write.csv(threshold3.qtl, file = paste0(output.files.prefix,".threshold3.qtl.csv"), row.names=FALSE, na="")

# Classify QTLs by LG and Peak Position
classified.qtl <- t.qtl[order(t.qtl$lg,t.qtl$pos.peak),]
classified.qtl$group <- with(classified.qtl,
                             paste0("chr",lg,"-mrk",marker))
write.csv(classified.qtl, file = paste0(output.files.prefix,".classified.qtl.csv"), row.names=FALSE, na="")
print("Done with QTL Analysis")
mpi.quit()
