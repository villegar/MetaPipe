x.normal.summary.mapping <- function(i){
  library(qtl)
  library(ggplot2)
  library(latex2exp)
  transformation.info <- normal.transformed.meansp$feature == features[i]
  transformation.info <- normal.transformed.meansp[transformation.info,c("transf","transf.value")][1,]
  record <- data.frame(
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
        new.record <- record[1,] # Create copy of record object
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
      if(k > 1){
        record <- rbind(record,new.record)
      }else{
        record <- new.record
      }
    }
    
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
  return(record)
}