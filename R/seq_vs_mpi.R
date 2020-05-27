FILES <- c("metabolomics.non.parametric.gen.csv","metabolomics.non.parametric.meansp.csv",
           "metabolomics.non.parametric.phe.csv","metabolomics.non.parametric.scanone.csv",
           "metabolomics.non.parametric.summary.mapping.csv",
           "metabolomics.normal.gen.csv","metabolomics.normal.meansp.csv",
           "metabolomics.normal.phe.csv","metabolomics.normal.scanone.csv",
           "metabolomics.summary.mapping.csv",
           "metabolomics.pareto.non.parametric.meansp.csv",
           "metabolomics.pareto.normal.meansp.csv",
           "metabolomics.transformed.meansp.csv")
MPIPATH <- "C:/Users/roberto.villegasdiaz/Downloads/metabolomics-pipeline/../QUIK/mpi-17/"
SEQPATH <- "C:/Users/roberto.villegasdiaz/Downloads/metabolomics-pipeline/../QUIK/190911/"

for(f in FILES){
  print(f)
  A <- read.csv(paste0(SEQPATH,f))
  B <- read.csv(paste0(MPIPATH,f))
  print(all.equal(A,B))
}

A <- read.csv(paste0(SEQPATH,"metabolomics.summary.mapping.csv"))
B <- read.csv(paste0(MPIPATH,"metabolomics.summary.mapping.csv"))
p5.lod <- cbind(A$p5.lod.thr,B$p5.lod.thr)
p5.lod <- p5.lod[!is.na(p5.lod[,1]),]
p5.lod <- p5.lod[!is.na(p5.lod[,2]),]
p5.lod[p5.lod[,1] != p5.lod[,2],]
cbind(as.character(A$qtl.ID),as.character(B$qtl.ID)) # 1 string mistmatch
# 395 : "M690T38:18@16.748043"   "M690T38:18@16.748044"
cbind(A$marker,B$marker) # 1 string mistmatch
cbind(A$p5.qtl,B$p5.qtl) # 5 mismatches

140
349
380
399
422