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
library(gplots)
library(grid)
library(gridExtra)
library(latex2exp)
library(MASS)
library(plyr)
library(psych)
library(R.devices)
library(tictoc)
library(tidyverse)
library(VIM)

# Parallel Libraries
library(foreach)
library(doParallel)
#library(Rmpi)
#library(snow)

# PCA Libraries
library(FactoMineR)
library(factoextra)

# QTL Analysis
library(qtl)

# Source util functions
# source("plots-util.R")
# source("transformations.R")

main <- function(){

# Obtain the number of available cores
cores <- parallel::detectCores()
CPUS <- cores[1] - 1

# closeAllConnections()

if(length(args) < 1){
  PERMUTATIONS <- 1000 # Number of permutations for QTL Analysis
  REPLACE_NA <- FALSE
  PARETO_SCALING <- FALSE
  OUT_PREFIX <- "metabolomics"
  PLOTS_DIR <- "metabolomics"
} else if(length(args) < 2){
  PERMUTATIONS <- as.numeric(args[1]) # Number of permutations for QTL Analysis
  REPLACE_NA <- FALSE
  PARETO_SCALING <- FALSE
  OUT_PREFIX <- "metabolomics"
  PLOTS_DIR <- "metabolomics"
} else if(length(args) < 3){
  PERMUTATIONS <- as.numeric(args[1]) # Number of permutations for QTL Analysis
  REPLACE_NA <- as.logical(args[2])
  PARETO_SCALING <- FALSE
  OUT_PREFIX <- "metabolomics"
  PLOTS_DIR <- "metabolomics"
} else if(length(args) < 4){
  PERMUTATIONS <- as.numeric(args[1]) # Number of permutations for QTL Analysis
  REPLACE_NA <- as.logical(args[2])
  PARETO_SCALING <- as.logical(args[3])
  OUT_PREFIX <- "metabolomics"
  PLOTS_DIR <- "metabolomics"
} else if(length(args) < 5){
  PERMUTATIONS <- as.numeric(args[1]) # Number of permutations for QTL Analysis
  REPLACE_NA <- as.logical(args[2])
  PARETO_SCALING <- as.logical(args[3])
  OUT_PREFIX <- args[4]
  PLOTS_DIR <- "metabolomics"
} else {
  PERMUTATIONS <- as.numeric(args[1]) # Number of permutations for QTL Analysis
  REPLACE_NA <- as.logical(args[2])
  PARETO_SCALING <- as.logical(args[3])
  OUT_PREFIX <- args[4]
  PLOTS_DIR <- args[5]
}

tictoc::tic.clearlog()
cat(paste0("CMD Parameters: (",PERMUTATIONS,",",REPLACE_NA,",",PARETO_SCALING,",",OUT_PREFIX,",",PLOTS_DIR,")"))

# Global parameters
excluded_columns <- c(1,2,3)
len_excluded_columns <- length(excluded_columns)
transf_vals <- c(2,exp(1))#,3,4,5,6,7,8,9,10
SEED <- 20190901 # Seed for QTL Analysis
LOD.THRESHOLD <- 3 # LOD threhold for QTL Analysis
prop_na <- 0.5 # Allows 50% of NAs per feature

# Environment configuration
dir.create(file.path(getwd(), PLOTS_DIR), showWarnings = FALSE) # Directory for plots

tictoc::tic("Total")
tictoc::tic("Loading and pre-processing")

# Load and Cleaning Data
raw_data <- MetaPipe::load_raw("sp.csv", excluded_columns)
raw_data_rows <- nrow(raw_data)

# Replacement of Missing Values
raw_data <- MetaPipe::replace_missing(raw_data, excluded_columns, OUT_PREFIX, prop_na, REPLACE_NA)

write.csv(raw_data, file = paste0(OUT_PREFIX,".all.raw_data.csv"), row.names=FALSE)
tictoc::toc(log = TRUE) # Loading and pre-processing

# generate.boxplots <- function(raw_data,ggplot_save){
#   print("Generating Boxplots")
#   cl <- makeCluster(CPUS, outfile=paste0('./info_parallel.log')) # Make cluster
#   registerDoParallel(cl)  # Register cluster
#   features <- colnames(raw_data)
#   AllPlots <- foreach(i=(len_excluded_columns + 1):ncol(raw_data), 
#                       .packages = c("ggplot2","latex2exp","R.devices")) %dopar% {
#                         myPlot <- ggplot(data=raw_data,aes(x=ID,y=raw_data[,i])) +
#                           geom_boxplot(aes(fill= "")) +
#                           theme(axis.text.x = element_text(angle = 60, hjust = 1))+ 
#                           labs(title=paste("Feature",features[i]), x='ID', y='')
#                         ggplot_save(myPlot,paste0("BOX_",(i - len_excluded_columns),"_",features[i]))
#                       }
#   stopCluster(cl1) # Stop cluster
#   print("Done with Boxplots")
# }
#generate.boxplots(raw_data,ggplot_save)

tictoc::tic("Normality Assessment")
features <- colnames(raw_data)
raw_data_normalised <- MetaPipe::assess_normality(raw_data, excluded_columns, CPUS, OUT_PREFIX, PLOTS_DIR, transf_vals)
# MetaPipe::assess_normality(raw_data, excluded_columns, cpus = 1, out_prefix = "metapipe", plots_dir = here::here("../testbed/"), transf_vals = c(2, exp(1)))
tictoc::toc(log = TRUE) # Normality Assessment
tictoc::tic("Transformed data post-processing")
#MetaPipe::assess_normality_postprocessing(raw_data, excluded_columns, raw_data_normalised, OUT_PREFIX, PARETO_SCALING)
#MetaPipe::assess_normality_stats(raw_data, excluded_columns, raw_data_normalised, OUT_PREFIX)
tictoc::toc(log = TRUE) # Transformed data post-processing

tictoc::tic("QTL analysis")
tictoc::tic("QTL analysis preprocessing")
# Prepocessing data for QTL Analysis
genetic_map <- read.csv("OriginalMap.csv")
colnames(genetic_map)[1] <- "ID"
genetic_map$ID <- as.character(genetic_map$ID)

## Normal features
colnames(raw_data_norm)[1] <- "ID"
raw_data_norm$GenoID <- with(raw_data_norm,
                             gsub(" ", "0", paste0(Generation, "_", sprintf("%3s", as.character(ID))))
                            )
raw_data_norm$ID <- raw_data_norm$GenoID
raw_data_norm$GenoID <- NULL

pheno_norm <- dplyr::inner_join(raw_data_norm, genetic_map, by = "ID")[, colnames(raw_data_norm)]
pheno_norm$Group <- NULL
pheno_norm$Generation <- NULL
geno_norm <- rbind(genetic_map[1:2, ], 
                   dplyr::inner_join(pheno_norm, genetic_map, by = "ID")[, colnames(genetic_map)]
                  )

## Non-parametric features
colnames(raw_data_non_par)[1] <- "ID"
raw_data_non_par$GenoID <- with(raw_data_non_par,
                                gsub(" ", "0", paste0(Generation, "_", sprintf("%3s", as.character(ID))))
                               )
raw_data_non_par$ID <- raw_data_non_par$GenoID
raw_data_non_par$GenoID <- NULL

pheno_non_par <- dplyr::inner_join(raw_data_non_par, genetic_map, by = "ID")[, colnames(raw_data_non_par)]
pheno_non_par$Group <- NULL
pheno_non_par$Generation <- NULL
geno_non_par <- rbind(genetic_map[1:2, ],
                      dplyr::inner_join(pheno_non_par, genetic_map, by = "ID")[, colnames(genetic_map)]
                     )

# Clean phenotypic data
empty_features_non_par <- sapply(pheno_non_par, function(x) all(is.na(x)) || all(is.infinite(x)))
empty_features_norm <- sapply(pheno_norm, function(x) all(is.na(x)) || all(is.infinite(x)))
#pheno_non_par.ncols <- ncol(pheno_non_par)
#pheno_norm.ncols <- ncol(pheno_norm)
pheno_non_par[empty_features_non_par] <- NULL
pheno_norm[empty_features_norm] <- NULL

if(any(empty_features_non_par)) {
  print(paste0("The following non-parametric features were removed (NAs):"))
  print(names(empty_features_non_par)[empty_features_non_par])
}

if(any(empty_features_norm)) {
  print(paste0("The following normal features were removed (NAs):"))
  print(names(empty_features_norm)[empty_features_norm])
}

# Write genotypic and phenotypic dataset
## Normal features
write.csv(geno_norm, file = paste0(OUT_PREFIX,"_geno_norm.csv"), row.names=FALSE)
write.csv(pheno_norm, file = paste0(OUT_PREFIX, "_pheno_norm.csv"), row.names=FALSE)
## Non-parametric features
write.csv(geno_non_par, file = paste0(OUT_PREFIX, "_geno_non_par.csv"), row.names=FALSE)
write.csv(pheno_non_par, file = paste0(OUT_PREFIX, "_pheno_non_par.csv"), row.names=FALSE)

tictoc::toc(log = TRUE) # QTL analysis preprocessing
tictoc::tic("Normal QTL analysis")
# QTL Analysis
x_norm <- qtl::read.cross("csvs", here::here(),
                          paste0(OUT_PREFIX, "_geno_norm.csv"),
                          paste0(OUT_PREFIX, "_pheno_norm.csv"))
features <- colnames(x_norm$pheno)
set.seed(SEED)
x_norm <- qtl::jittermap(x_norm)
x_norm <- qtl::calc.genoprob(x_norm, step = 1, error.prob = 0.001)
num_indv_phend_n <- summary(x_norm)[[2]]
print("Starting with QTL Analysis")
print("Starting with Normal QTL Analysis")

# Obtain LOD scores for all features and markers
# cl <- parallel::makeCluster(ceiling(CPUS*1), outfile=paste0('./info_parallel_QTL.log'))
# doParallel::registerDoParallel(cl)
# x_norm_scone <- foreach(i=2:ncol(x_norm$pheno),
#                      .combine = cbind) %dopar% {
#                        
#                        # Run single scan
#                        normal.scanone <- qtl::scanone(x_norm, pheno.col = i,  model = "normal", method = "hk")
#                        if(i == 2){
#                          record <- data.frame(
#                            chr = normal.scanone$chr,
#                            pos = normal.scanone$pos,
#                            lod = normal.scanone$lod,
#                            row.names = rownames(normal.scanone)
#                          )
#                          colnames(record)[3] <- features[i]
#                        }
#                        else{
#                          record <- data.frame(data = normal.scanone$lod)
#                          colnames(record) <- features[i]
#                        }
#                        record
#                      }
# parallel::stopCluster(cl) # Stop cluster
x_norm_scone <- MetaPipe::qtl_scone(x_norm, CPUS)
write.csv(x_norm_scone, file = paste0(OUT_PREFIX, "_normal_scanone.csv"))

cl <- parallel::makeCluster(ceiling(CPUS*0.5), outfile=paste0('./info_parallel_QTL.log'))
doParallel::registerDoParallel(cl)
x_norm_sum_map <- foreach(i=2:ncol(x_norm$pheno),
                         .combine = rbind) %dopar% {
                           transformation.info <- raw_data_transformed_norm$feature == features[i]
                           transformation.info <- raw_data_transformed_norm[transformation.info,c("transf","transf.value")][1,]
                           
                           record <- data.frame(
                             ID = i - 1,
                             qtl.ID = NA,
                             trait = features[i],
                             ind = num_indv_phend_n,
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
                               marker.info <- qtl::find.markerpos(cross, 
                                                                  qtl::find.marker(cross, chr = chr, pos = pos))
                               new.marker <- rownames(marker.info)
                               new.pos <- marker.info$pos
                             }
                             return(c(new.marker,as.character(new.pos)))
                           }
                           
                           # Run single scan
                           normal.scanone <-  qtl::scanone(x_norm, pheno.col = i,  model = "normal", method = "hk")
                           summary.normal.scanone <- summary(normal.scanone, threshold = LOD.THRESHOLD)
                           lod.count <- nrow(summary.normal.scanone)
                           if(!is.null(lod.count) && lod.count > 0) {
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
                               marker.info <- transform.pseudomarker(x_norm,marker,new.record$lg,new.record$pos.peak)
                               new.record$marker <- marker.info[1]
                               new.record$pos.peak <- as.numeric(marker.info[2])
                               
                               if(!is.na(new.record$lg)){
                                 new.record$qtl.ID <- with(new.record, sprintf("%s:%s@%f",features[i],lg,pos.peak))
                               }
                               
                               p95.bayesian <- qtl::bayesint(normal.scanone, chr = new.record$lg ,expandtomarkers = TRUE, prob = 0.95)
                               p95.bayesian <- unique(p95.bayesian)
                               #p95.bayesian <- summary(normal.scanone,  perms=normal.scanone.per, alpha=0.5, pvalues=TRUE)
                               low.bound <- 1#p95.bayesian$pos == min(p95.bayesian$pos)
                               upper.bound <- p95.bayesian$pos == max(p95.bayesian$pos)
                               
                               p95.bayesian$marker <- NA # Add new column for markers, prevent duplicated row names
                               # Verify if the Bayesian interval QTLs have pseudomarkers
                               for(l in 1:nrow(p95.bayesian)){
                                 marker <- rownames(p95.bayesian)[l]
                                 marker.info <- transform.pseudomarker(x_norm,marker,p95.bayesian[l,"chr"],p95.bayesian[l,"pos"])
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
                               normal.scanone.per <- qtl::scanone(x_norm, pheno.col = i, model = "normal", method = "hk", n.perm = PERMUTATIONS)
                               p5 <- summary(normal.scanone.per)[[1]]  #  5% percent
                               p10 <- summary(normal.scanone.per)[[2]] # 10% percent
                               
                               
                               lod.plot <- MetaPipe::save_plot(plot(normal.scanone, ylab="LOD Score") + 
                                                      abline(h=p5, lwd=2, lty="solid", col="red") +
                                                      abline(h=p10, lwd=2, lty="solid", col="red"),
                                                    paste0(PLOTS_DIR,"/LOD-",features[i]), width = 18)
                               
                               record[,]$p5.lod.thr <- p5
                               record[,]$p10.lod.thr <- p10
                               
                               p5.index <- record$lod.peak >= p5
                               p10.index <- record$lod.peak >= p10
                               if(!is.na(p5.index) && any(p5.index)){ record[p5.index,]$p5.qtl <- TRUE }
                               if(!is.na(p10.index)&& any(p10.index)){ record[p10.index,]$p10.qtl <- TRUE }
                               

                               chr <- as.numeric(summary.normal.scanone$chr)
                               pos <- as.numeric(summary.normal.scanone$pos)
                               qtl_s <- qtl::makeqtl(x_norm, chr, pos, what=c("prob"))
                               
                               for(m in 1:length(chr)){
                                 #qtl_s <- makeqtl(x_norm, chr[m], pos[m], what=c("prob"))
                                 #f <- as.formula(paste0("y~",paste0("Q",seq(1:nrow(summary.normal.scanone)), collapse = " + ")))
                                 f <- as.formula(paste0("y~",paste0("Q",m, collapse = " + ")))
                                 fitqtl <- qtl::fitqtl(x_norm, pheno.col = i, qtl_s, formula = f , get.ests = TRUE, model = "normal", method="hk")
                                 summary.fitqtl <- summary(fitqtl)
                                 
                                 if(length(summary.fitqtl)){
                                   p.var <- as.numeric(summary.fitqtl[[1]][1,"%var"])
                                   p.value.f <- as.numeric(summary.fitqtl[[1]][,"Pvalue(F)"])[1]
                                   estimates <- as.numeric(summary.fitqtl$ests[,"est"])[-1]
                                   record[m,]$pvar <- p.var
                                   record[m,]$p.val <- p.value.f
                                   record[m,]$est.add <- estimates[1]
                                   record[m,]$est.dom <- estimates[2]
                                   #for(l in 1:length(estimates)){
                                   #  offset <- 2*(l-1)
                                   #  record[l,]$est.add <- estimates[offset + 1]
                                   #  record[l,]$est.dom <- estimates[offset + 2]
                                   #}
                                 }
                               }
                               
                               # No needed for this data set
                               #refinqtl <- refineqtl(x_norm, qtl = qtl_s, pheno.col = i, formula = f, verbose = FALSE, model = "normal", method="hk")
                               #refinqtl
                               
                               #fitqtl <- fitqtl(x_norm, pheno.col = i, refinqtl, formula = f, get.ests = TRUE, model = "normal", method="hk")
                               #summary(fitqtl)
                               
                               
                               ## find additional QTLs
                               #out.aq <- addqtl(x_norm, qtl = refinqtl, pheno.col = i, formula = f, method="hk")
                               #max(out.aq)
                           }
                           record
                         }
parallel::stopCluster(cl) # Stop cluster
#tmp <- x_norm_sum_map
#tmp$qtl[!is.na(tmp$qtl)] <- 1:length(tmp$qtl[!is.na(tmp$qtl)])
#x_norm_sum_map <- tmp
write.csv(x_norm_sum_map, file = paste0(OUT_PREFIX, "_norm_summ_map.csv"), row.names = FALSE, na = "")

tictoc::toc(log = TRUE) # Normal QTL analysis
tictoc::tic("Non-parametric QTL analysis")
# Non-parametric QTL
x_non_par <- qtl::read.cross("csvs", here::here(),
                             paste0(OUT_PREFIX, "_geno_non_par.csv"),
                             paste0(OUT_PREFIX, "_pheno_non_par.csv"))
features_np <- colnames(x_non_par$pheno)
set.seed(SEED)
x_non_par <- qtl::jittermap(x_non_par)
x_non_par <- qtl::calc.genoprob(x_non_par, step = 1, error.prob = 0.001)
num_indv_phend_np <- summary(x_non_par)[[2]]
print("Starting with Non-Parametric QTL Analysis")

# Obtain LOD scores for all features and markers
cl <- parallel::makeCluster(ceiling(CPUS*1), outfile=paste0('./info_parallel_QTL.log'))
doParallel::registerDoParallel(cl)
x_non_par_scone <- foreach(i=2:ncol(x_non_par$pheno),
                     .combine = cbind) %dopar% {
                       
                       # Run single scan
                       non.parametric.scanone <- qtl::scanone(x_non_par, pheno.col = i,  model = "np")
                       if(i == 2){
                         record <- data.frame(
                           chr = non.parametric.scanone$chr,
                           pos = non.parametric.scanone$pos,
                           lod = non.parametric.scanone$lod,
                           row.names = rownames(non.parametric.scanone)
                         )
                         colnames(record)[3] <- features_np[i]
                       }
                       else{
                         record <- data.frame(data = non.parametric.scanone$lod)
                         colnames(record) <- features_np[i]
                       }
                       record
                     }
parallel::stopCluster(cl) # Stop cluster
write.csv(x_non_par_scone, file = paste0(OUT_PREFIX,"_non_par_scanone.csv"))

cl <- parallel::makeCluster(ceiling(CPUS*0.5), outfile=paste0('./info_parallel_QTL.log'))
doParallel::registerDoParallel(cl)
x_non_par_sum_map <- foreach(i=2:ncol(x_non_par$pheno),
                             .combine = rbind) %dopar% {
                               #transformation.info <- raw_data_transformed_non_par$feature == features_np[i]
                               #transformation.info <- raw_data_transformed_non_par[transformation.info,c("transf","transf.value")][1,]
                               
                               record <- data.frame(
                                 ID = i - 1,
                                 qtl.ID = NA,
                                 trait = features_np[i],
                                 ind = num_indv_phend_np,
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
                                   marker.info <- qtl::find.markerpos(cross, 
                                                                      qtl::find.marker(cross, chr = chr, pos = pos))
                                   new.marker <- rownames(marker.info)
                                   new.pos <- marker.info$pos
                                 }
                                 return(c(new.marker,as.character(new.pos)))
                               }
                               
                               # Run single scan
                               non.parametric.scanone <- qtl::scanone(x_non_par, pheno.col = i,  model = "np")
                               summary.non.parametric.scanone <- summary(non.parametric.scanone, threshold = LOD.THRESHOLD)
                               lod.count <- nrow(summary.non.parametric.scanone)
                               if(!is.null(lod.count) && lod.count > 0) {
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
                                   marker.info <- transform.pseudomarker(x_non_par,marker,new.record$lg,new.record$pos.peak)
                                   new.record$marker <- marker.info[1]
                                   new.record$pos.peak <- as.numeric(marker.info[2])
                                   
                                   if(!is.na(new.record$lg)){
                                     new.record$qtl.ID <- with(new.record, sprintf("%s:%s@%f",features_np[i],lg,pos.peak))
                                   }
                                   
                                   p95.bayesian <- qtl::bayesint(non.parametric.scanone, chr = new.record$lg, expandtomarkers = TRUE, prob = 0.95)
                                   p95.bayesian <- unique(p95.bayesian)
                                   #p95.bayesian <- summary(non.parametric.scanone,  perms=non.parametric.scanone.per, alpha=0.5, pvalues=TRUE)
                                   low.bound <- 1#p95.bayesian$pos == min(p95.bayesian$pos)
                                   upper.bound <- p95.bayesian$pos == max(p95.bayesian$pos)
                                   
                                   p95.bayesian$marker <- NA # Add new column for markers, prevent duplicated row names
                                   # Verify if the Bayesian interval QTLs have pseudomarkers
                                   for(l in 1:nrow(p95.bayesian)){
                                     marker <- rownames(p95.bayesian)[l]
                                     marker.info <- transform.pseudomarker(x_non_par,marker,p95.bayesian[l,"chr"],p95.bayesian[l,"pos"])
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
                                 
                                 non.parametric.scanone.per <- qtl::scanone(x_non_par, pheno.col = i, model = "np", n.perm = PERMUTATIONS)
                                 p5 <- summary(non.parametric.scanone.per)[[1]]  #  5% percent
                                 p10 <- summary(non.parametric.scanone.per)[[2]] # 10% percent
                                 
                                 
                                 lod.plot <- MetaPipe::save_plot(plot(non.parametric.scanone, ylab="LOD Score") + 
                                                        abline(h=p5, lwd=2, lty="solid", col="red") +
                                                        abline(h=p10, lwd=2, lty="solid", col="red"),
                                                      paste0(PLOTS_DIR,"/LOD-NP-", features_np[i]), width = 18)
                                 
                                 record[,]$p5.lod.thr <- p5
                                 record[,]$p10.lod.thr <- p10
                                 
                                 p5.index <- record$lod.peak >= p5
                                 p10.index <- record$lod.peak >= p10
                                 if(!is.na(p5.index) && any(p5.index)){ record[p5.index,]$p5.qtl <- TRUE }
                                 if(!is.na(p10.index)&& any(p10.index)){ record[p10.index,]$p10.qtl <- TRUE }
                               }
                               record
                             }
parallel::stopCluster(cl) # Stop cluster

write.csv(x_non_par_sum_map, file = paste0(OUT_PREFIX, "_non_par_summ_map.csv"), row.names=FALSE, na="")

tictoc::toc(log = TRUE) # Non-parametric QTL analysis
tictoc::tic("QTL analysis postprocessing")
x_norm <- qtl::read.cross("csvs", here::here(),
                          paste0(OUT_PREFIX, "_geno_norm.csv"),
                          paste0(OUT_PREFIX, "_pheno_norm.csv"))
features <- colnames(x_norm$pheno)
set.seed(SEED)
x_norm <- qtl::jittermap(x_norm)
x_norm <- qtl::calc.genoprob(x_norm, step=1, error.prob=0.001)
num_indv_phend_n <- summary(x_norm)[[2]]

# Non-parametric QTL
x_non_par <- qtl::read.cross("csvs", here::here(),
                             paste0(OUT_PREFIX, "_geno_non_par.csv"),
                             paste0(OUT_PREFIX, "_pheno_non_par.csv"))
features_np <- colnames(x_non_par$pheno)
set.seed(SEED)
x_non_par <- qtl::jittermap(x_non_par)
x_non_par <- qtl::calc.genoprob(x_non_par, step=1, error.prob=0.001)
num_indv_phend_np <- summary(x_non_par)[[2]]

# Generate effect plots
x_non_par_sim <- qtl::sim.geno(x_non_par)
x_norm_sim <- qtl::sim.geno(x_norm)

true_qtl <- plyr::rbind.fill(x_norm_sum_map[x_norm_sum_map$p5.qtl,],
                             x_non_par_sum_map[x_non_par_sum_map$p5.qtl,])
thrsh3_qtl <- plyr::rbind.fill(x_norm_sum_map[x_norm_sum_map$lod.peak > LOD.THRESHOLD,],
                               x_non_par_sum_map[x_non_par_sum_map$lod.peak > LOD.THRESHOLD,])
true_qtl_features <- as.character(true_qtl$trait)
true_qtl_markers <- as.character(true_qtl$marker)

cl <- parallel::makeCluster(ceiling(CPUS), outfile=paste0('./info_parallel_QTL.log'))
doParallel::registerDoParallel(cl)
effect_plots <- foreach(i=1:nrow(true_qtl),
                        .packages = c("latex2exp","qtl","R.devices")) %dopar% {
                          if(true_qtl[i,]$method == "normal-scanone"){
                            if(true_qtl[i,]$transf == "log"){
                              ylab <- paste0("$\\log_{",true_qtl[i,]$transf.val,"}(",true_qtl_features[i],")$")
                            } else if(true_qtl[i,]$transf == "root"){
                              ylab <- paste0("$\\sqrt[",true_qtl[i,]$transf.val,"]{",true_qtl_features[i],"}$")
                            } else if(true_qtl[i,]$transf == "power"){
                              ylab <- paste0("$(",true_qtl_features[i],")^",true_qtl[i,]$transf.val,"$")
                            } else {
                              ylab <- true_qtl_features[i]
                            }
                            effect_plots <- MetaPipe::save_plot(qtl::effectplot(x_norm_sim, pheno.col = true_qtl_features[i], 
                                                                      mname1 = true_qtl_markers[i], main = NULL, ylab = latex2exp::TeX(ylab)),
                                                      paste0(PLOTS_DIR,"/EFF-",true_qtl_features[i],"-",true_qtl_markers[i]))
                          } else {
                            ylab <- true_qtl_features[i]
                            effect_plots <- MetaPipe::save_plot(qtl::effectplot(x_non_par_sim, pheno.col = as.character(true_qtl_features[i]), 
                                                                      mname1 = true_qtl_markers[i], main = NULL, ylab = latex2exp::TeX(ylab)),
                                                      paste0(PLOTS_DIR,"/EFF-NP-",true_qtl_features[i],"-",true_qtl_markers[i]))
                          }
                        }
parallel::stopCluster(cl) # Stop cluster

write.csv(true_qtl, file = paste0(OUT_PREFIX,".true.qtl.csv"), row.names=FALSE, na="")
write.csv(thrsh3_qtl, file = paste0(OUT_PREFIX,".threshold3.qtl.csv"), row.names=FALSE, na="")

# Classify QTLs by LG and Peak Position
class_qtl <- true_qtl[order(true_qtl$lg,true_qtl$pos.peak),]
class_qtl$group <- with(class_qtl,
                        paste0("chr",lg,"-mrk",marker))
write.csv(class_qtl, file = paste0(OUT_PREFIX,".classified.qtl.csv"), row.names=FALSE, na="")
print("Done with QTL Analysis")
tictoc::toc(log = TRUE) # QTL analysis postprocessing
tictoc::toc(log = TRUE) # QTL analysis


# For both PCA and LDA the data must have no NAs and must be scaled
raw_data <- read.csv(paste0(OUT_PREFIX,".all.raw_data.csv"))
if(!REPLACE_NA){
  NA2halfmin <- function(x) suppressWarnings(replace(x, is.na(x), (min(x, na.rm = TRUE)/2)))
  raw_data[,-excluded_columns] <- lapply(raw_data[,-excluded_columns], NA2halfmin)
}
raw_data_transformed <- raw_data # No scaling
#if(!PARETO_SCALING){ # Apply Pareto Scaling
#  raw_data_transformed <- cbind(raw_data[,excluded_columns],paretoscale(raw_data[,-excluded_columns]))
#}

tictoc::tic("PCAnalysis")
# PCAnalysis with mean (used no missing data) 
##OUT_PREFIX <- "S1-metabolomics"
##raw_data_transformed <- read.csv(paste0(OUT_PREFIX,".all.raw_data.csv"))
##raw_data_transformed$X <- NULL
##raw_data_transformed <- raw_data_transformed[order(as.character(raw_data_transformed$ID)),]
##raw_data_transformed$Group <- NULL
res.pca <- PCA(raw_data_transformed[,-excluded_columns],  graph = FALSE, scale.unit = TRUE)
##fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
##res.pca$eig
# Biplot with top 10 features 
save_plotTIFF(fviz_pca_biplot(res.pca, col.var="contrib",
                         gradient.cols = c("green","red","blue"),#"#00AFBB" "#E7B800", "#FC4E07"),
                         select.var = list(contrib = 10),
                         label="var",addEllipses=TRUE, ellipse.level=0.95, repel = TRUE  # Avoid text overlapping
),
paste0(PLOTS_DIR,"/PCA-biplot.top10"),12,6)
tictoc::toc(log = TRUE) # PCAnalysis
tictoc::tic("LDAnalysis")
# LDAnalysis
## Create an "unknown" group name for missing data
raw_data_transformed$Group <- as.character(raw_data_transformed$Group)
raw_data_transformed$Group[is.na(raw_data_transformed$Group)] <- "Unknown"

## Calculate mean by color
raw_data_transformed.diff.by.color <- 
  data.frame(t(aggregate(raw_data_transformed[,-excluded_columns], list(raw_data_transformed$Group), mean))[-1,])
raw_data_transformed.diff.by.color$X1 <- as.numeric(as.character(raw_data_transformed.diff.by.color$X1))
raw_data_transformed.diff.by.color$X2 <- as.numeric(as.character(raw_data_transformed.diff.by.color$X2))
colnames(raw_data_transformed.diff.by.color) <- c("black.mean","white.mean","unknown.mean")
raw_data_transformed.diff.by.color$mean.diff <- with(raw_data_transformed.diff.by.color, black.mean-white.mean)
raw_data_transformed.diff.by.color <- raw_data_transformed.diff.by.color[order(raw_data_transformed.diff.by.color$black.mean, decreasing = T),]
top.100.black <- rownames(raw_data_transformed.diff.by.color)[1:100]
raw_data_transformed.diff.by.color <- raw_data_transformed.diff.by.color[order(raw_data_transformed.diff.by.color$white.mean, decreasing = T),]
top.100.white <- rownames(raw_data_transformed.diff.by.color)[1:100]
## Whole dataset and Top 200 features LDA

top.200 <- unique(c(top.100.black,top.100.white))
#top.200 <- rownames(raw_data_transformed.diff.by.color)[1:200] # There's something funny after 75
colored.raw_data_transformed.full <- cbind(raw_data_transformed$Group,raw_data_transformed[,-excluded_columns])
colored.raw_data_transformed.top200 <- cbind(raw_data_transformed$Group,raw_data_transformed[,top.200])
colnames(colored.raw_data_transformed.full)[1] <- "FruitColor"
colnames(colored.raw_data_transformed.top200)[1] <- "FruitColor"
fit.full <- lda(FruitColor ~ ., data = colored.raw_data_transformed.full)
fit.top200 <- lda(FruitColor ~ ., data = colored.raw_data_transformed.top200)

lda.data.full <- cbind(colored.raw_data_transformed.full, predict(fit.full)$x)
lda.data.top200 <- cbind(colored.raw_data_transformed.top200, predict(fit.top200)$x)
save_plotTIFF(ggplot(lda.data.full, aes(LD1,LD2)) +
  geom_point(aes(color = FruitColor)) +
  stat_ellipse(aes(x=LD1, y=LD2, fill = FruitColor), alpha = 0.2, geom = "polygon"),
  paste0(PLOTS_DIR,"/LDA-full-dataset"))
  
save_plotTIFF(ggplot(lda.data.top200, aes(LD1,LD2)) +
  geom_point(aes(color = FruitColor)) +
  stat_ellipse(aes(x=LD1, y=LD2, fill = FruitColor), alpha = 0.2, geom = "polygon"),
  paste0(PLOTS_DIR,"/LDA-top200-dataset.h7"), height = 7)
tictoc::toc(log = TRUE) # LDAnalysis

tictoc::tic("Heatmap for true QTLs")
# Heatmap
x_norm.lod.scores <- read.csv(paste0(OUT_PREFIX,".normal.scanone.csv"))
x.non.parametric.lod.scores <- read.csv(paste0(OUT_PREFIX,".non.parametric.scanone.csv"))
true.qtl <- read.csv(paste0(OUT_PREFIX,".true.qtl.csv"))
true.qtl.features <- unique(as.character(true.qtl$trait))
x_norm.lod.scores.true.qtl <- x_norm.lod.scores[, which(names(x_norm.lod.scores) %in% true.qtl.features)]
x.non.parametric.lod.scores.true.qtl <- x.non.parametric.lod.scores[, which(names(x.non.parametric.lod.scores) %in% true.qtl.features)]
lod.scores <- cbind(x_norm.lod.scores.true.qtl,x.non.parametric.lod.scores.true.qtl)
lod.scores <- matrix(as.numeric(unlist(lod.scores)), nrow=nrow(lod.scores))
rownames(lod.scores) <- x_norm.lod.scores$X
colnames(lod.scores) <- true.qtl.features
obs.by.chr <- table(x_norm.lod.scores$chr)
colnams.chr <- rep(NA,length(x_norm.lod.scores$X))
k <- 1
for(i in 1:length(obs.by.chr)){
  colnams.chr[k] <- paste0("chr ",i)
  k <- k + obs.by.chr[i]
}
#rownames(lod.scores) <- colnams.chr

# Heatmap without dendrogram
save_plot(heatmap.2(lod.scores, Rowv = FALSE, Colv = FALSE, 
                   scale = "none",
                   margins = c(6, 1),
                   trace = "none", 
                   tracecol = "black",
                   symkey = FALSE, 
                   symbreaks = FALSE, 
                   dendrogram = "none",
                   density.info = "histogram", 
                   denscol = "black",
                   key.par=list(mar=c(3.5,0,3,0)),
                   col = rev(heat.colors(n = 10)),
                   cexRow = 0.8,
                   labRow = colnams.chr,
                   lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2.5, 5), lwid=c(1, 10, 1))
         ,paste0(PLOTS_DIR,"/HEAT-without-dendrogram", 16, 16))

# Heatmap with dendrogram and key at the top-left corner
save_plotTIFF(heatmap.2(lod.scores, Rowv = FALSE, Colv = TRUE, 
                       scale = "none",
                       margins = c(6, 6),
                       trace = "none", 
                       tracecol = "black",
                       symkey = FALSE, 
                       symbreaks = FALSE, 
                       dendrogram = "column",
                       density.info = "histogram", 
                       denscol = "black",
                       col = rev(heat.colors(n = 100)),
                       cexRow = 0.8,
                       cexCol = 0.8,
                       labRow = colnams.chr,
                       lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2.5, 5), lwid=c(1, 10, 1),
                       key.par=list(mar=c(3.5,3,3,0)))
             ,paste0(PLOTS_DIR,"/HEAT-with-dendrogram-all.100cov2"))

# Heatmap with dendrogram and key at the bottom
save_plotTIFF(heatmap.2(lod.scores, Rowv = FALSE, Colv = TRUE, 
                       scale = "none",
                       #margins = c(6, 6),
                       trace = "none", 
                       tracecol = "black",
                       symkey = FALSE, 
                       symbreaks = FALSE, 
                       dendrogram = "column",
                       density.info = "histogram", 
                       denscol = "black",
                       col = rev(heat.colors(n = 100)),
                       cexRow = 0.8,
                       cexCol = 0.8,
                       labRow = colnams.chr,
                       lmat=rbind(c(0,3),c(2,1),c(0,4)), lhei=c(2,6,2), lwid=c(0.3, 7),
                       key.par=list(mar=c(3.5,1.5,2.5,5)))
             ,paste0(PLOTS_DIR,"/HEAT-with-dendrogram-all-bottom-key.100co"))
tictoc::toc(log = TRUE) # Heatmap for true QTLs
# closeAllConnections()
tictoc::toc(log = TRUE) # Total
log.txt <- tictoc::tic.log(format = TRUE)
write(unlist(log.txt), paste0(OUT_PREFIX,".log.times.p",PERMUTATIONS,".txt"))
}