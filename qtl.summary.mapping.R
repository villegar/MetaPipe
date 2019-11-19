args = commandArgs(trailingOnly=TRUE) # Command line arguments
# General Libraries
library(Amelia, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(gplots)
library(grid, quietly = TRUE)
library(gridExtra, quietly = TRUE)
library(latex2exp, quietly = TRUE)
library(MASS, quietly = TRUE)
library(plyr, quietly = TRUE)
library(psych, quietly = TRUE)
library(R.devices, quietly = TRUE)
library(tictoc, quietly = TRUE)
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
source("transformations.R")

# Resources
cores <- detectCores()
CPUS <- cores[1] - 1

if(length(args) < 1){
  PERMUTATIONS <- 1000 # Number of permutations for QTL Analysis
  REPLACE.NA <- FALSE
  PARETO.SCALING <- FALSE
  OUT.PREFIX <- "metabolomics"
  PLOTS.DIR <- "metabolomics"
} else if(length(args) < 2){
  PERMUTATIONS <- as.numeric(args[1]) # Number of permutations for QTL Analysis
  REPLACE.NA <- FALSE
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

tic.clearlog()
# Global parameters
SEED <- 20190901 # Seed for QTL Analysis
LOD.THRESHOLD <- 3 # LOD threhold for QTL Analysis

tic("Normal QTL analysis: Summary mapping")
# QTL Summary Mapping
# Load MPI libraries
library(Rmpi)
library(doMPI)

# In case R exits unexpectedly, have it automatically clean up
# resources taken up by Rmpi (slaves, memory, etc...)
.Last <- function(){
  if (is.loaded("mpi_initialize")){
    if (mpi.comm.size(1) > 0){
      print("Please use mpi.close.Rslaves() to close workers.")
      mpi.close.Rslaves()
    }
    print("Please use mpi.quit() to quit R")
    .Call("mpi_finalize")
  }
}

# Obtain LOD scores for all features and markers
cl <- startMPIcluster()
registerDoMPI(cl)

# QTL Analysis
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

individuals.phenotyped.np <- summary(x.non.parametric)[[2]]

normal.transformed.meansp <- read.csv(paste0(OUT.PREFIX,".normal.transformations.summary.csv"))

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
                                                             paste0(PLOTS.DIR,"/LOD-",features[i]), width = 18)
                                        
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
toc(log = TRUE) # Normal QTL analysis: Summary mapping
tic("Non-parametric QTL analysis: Summary mapping")
# Non-parametric QTL
print("Starting with Non-Parametric QTL Analysis")

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
                                                                     paste0(PLOTS.DIR,"/LOD-NP-",features.np[i]), width = 18)
                                                
                                                record[,]$p5.lod.thr <- p5
                                                record[,]$p10.lod.thr <- p10
                                                
                                                p5.index <- record$lod.peak >= p5
                                                p10.index <- record$lod.peak >= p10
                                                if(!is.na(p5.index) && any(p5.index)){ record[p5.index,]$p5.qtl <- TRUE }
                                                if(!is.na(p10.index)&& any(p10.index)){ record[p10.index,]$p10.qtl <- TRUE }
                                                
                                              }
                                              record
                                            }
toc(log = TRUE) # Non-normal QTL analysis: Summary mapping
tic("Effect plots and QTL analysis postprocessing")
# Generate effect plots
x2.non.parametric <- sim.geno(x.non.parametric)
x2.normal <- sim.geno(x.normal)
#effectplot(x2, pheno.col = "M155T28", mname1 = "gbs_13_305342", main = NULL)
t.qtl <- rbind.fill(x.normal.summary.mapping[x.normal.summary.mapping$p5.qtl,],
                    x.non.parametric.summary.mapping[x.non.parametric.summary.mapping$p5.qtl,])
threshold3.qtl <- rbind.fill(x.normal.summary.mapping[x.normal.summary.mapping$lod.peak > LOD.THRESHOLD,],
                             x.non.parametric.summary.mapping[x.non.parametric.summary.mapping$lod.peak > LOD.THRESHOLD,])
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
                                                    paste0(PLOTS.DIR,"/EFF-",features.t.qtl[i],"-",markers.t.qtl[i]))
                          } else {
                            ylab <- features.t.qtl[i]
                            effect.plot <- savePlot(effectplot(x2.non.parametric, pheno.col = as.character(features.t.qtl[i]), 
                                                               mname1 = markers.t.qtl[i], main = NULL, ylab = TeX(ylab)),
                                                    paste0(PLOTS.DIR,"/EFF-NP-",features.t.qtl[i],"-",markers.t.qtl[i]))
                          }
                        }

closeCluster(cl) # Stop cluster
#stopCluster.MPIcluster(c1)

#write.csv(x.normal.scanone, file = paste0(OUT.PREFIX,".normal.scanone.csv"))
write.csv(x.normal.summary.mapping, file = paste0(OUT.PREFIX,".normal.summary.mapping.csv"), row.names=FALSE, na="")
#write.csv(x.non.parametric.scanone, file = paste0(OUT.PREFIX,".non.parametric.scanone.csv"))
write.csv(x.non.parametric.summary.mapping, file = paste0(OUT.PREFIX,".non.parametric.summary.mapping.csv"), row.names=FALSE, na="")
write.csv(t.qtl, file = paste0(OUT.PREFIX,".true.qtl.csv"), row.names=FALSE, na="")
write.csv(threshold3.qtl, file = paste0(OUT.PREFIX,".threshold3.qtl.csv"), row.names=FALSE, na="")

# Classify QTLs by LG and Peak Position
classified.qtl <- t.qtl[order(t.qtl$lg,t.qtl$pos.peak),]
classified.qtl$group <- with(classified.qtl,
                             paste0("chr",lg,"-mrk",marker))
write.csv(classified.qtl, file = paste0(OUT.PREFIX,".classified.qtl.csv"), row.names=FALSE, na="")
toc(log = TRUE) # Effect plots and QTL analysis postprocessing
print("Done with QTL Analysis")
toc(log = TRUE) # QTL analysis
