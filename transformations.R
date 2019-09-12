# Log Transformation and Pareto scaling
## Pareto scaling function
#This function is adapted from Stephen C. Grace and Dane A. Hudson 
paretoscale2=function(z){
  colmean=apply(z,2,mean)
  colsd=apply(z,2,sd)
  colsqrtsd=sqrt(colsd)
  rv=sweep(z,2,colmean,"-")
  rv=sweep(rv,2,colsqrtsd, "/")
  return(rv)
}

paretoscale <- function(z){
  return(z)
}

## Log transformation
log.transformation <- function(shapiro, data, feature){
  logBases <- c(2,exp(1),3,4,5,6,7,8,9,10)
  
  record <- data.frame(
    index = i,
    feature = features[i],
    values = meansp[,i],
    flag = "Non-normal",
    transf = "",
    transf.value = 0
  )
  
  for(base in logBases){
    pvalue <- shapiro.test(log(data,base))[[2]]
    if(pvalue >= 0.05){
      transformed <- log(data,base)
      if(base == exp(1))
        base <- "e"
      #print(paste0("After log_",base," transformation ",pvalue," - originally ",shapiro))
      xlab <- paste0("$\\log_{",base,"}(",feature,")$")
      transformation <- paste0("LOG_",base)
      name.prefix <- paste0("plots/HIST_",(i-3),"_",transformation)
      generate.overlay.histogram(data,transformed,feature,name.prefix,xlab)
      record$values <- transformed
      record$transf <- "log"
      record$transf.value <- base
      return(record)
    }
  }
  return(root.transformation(shapiro, data, feature))
}

## Power transformation
power.transformation <- function(shapiro, data, feature){
  powers <- c(2,exp(1),3,4,5,6,7,8,9,10)
  
  record <- data.frame(
    index = i,
    feature = features[i],
    values = meansp[,i],
    flag = "Non-normal",
    transf = "",
    transf.value = 0
  )
  
  for(p in powers){
    pvalue <- shapiro.test(data^p)[[2]]
    if(pvalue >= 0.05){
      transformed <- data^p
      if(p == exp(1))
        p <- "e"
      #print(paste0("After X^",p," transformation ",pvalue," - originally ",shapiro))
      #is.pow <- TRUE
      xlab <- paste0("$(",feature,")^",p,"$")
      transformation <- paste0("POW_",p)
      name.prefix <- paste0("plots/HIST_",(i-3),"_",transformation)
      generate.overlay.histogram(data,transformed,feature,name.prefix,xlab)
      record$values <- transformed
      record$transf <- "power"
      record$transf.value <- p
      return(record)
    }
  }
  return(data.frame())
}

## Root transformation
root.transformation <- function(shapiro, data, feature){
  roots <- c(2,exp(1),3,4,5,6,7,8,9,10)
  
  record <- data.frame(
    index = i,
    feature = features[i],
    values = meansp[,i],
    flag = "Non-normal",
    transf = "",
    transf.value = 0
  )
  
  for(r in roots){
    pvalue <- shapiro.test(data^(1/r))[[2]]
    if(pvalue >= 0.05){
      transformed <- data^(1/r)
      if(r == exp(1))
        r <- "e"
      #print(paste0("After X^(1/",r,") transformation ",pvalue," - originally ",shapiro))
      xlab <- paste0("$\\sqrt[",r,"]{",feature,"}$")
      transformation <- paste0("ROOT_",r)
      name.prefix <- paste0("plots/HIST_",(i-3),"_",transformation)
      generate.overlay.histogram(data,transformed,feature,name.prefix,xlab)
      record$values <- transformed
      record$transf <- "root"
      record$transf.value <- r
      return(record)
    }
  }
  return(power.transformation(shapiro, data, feature, record))
}

## 
transform.data <- function(shapiro, data, feature, index, 
                           offset = 3, plots.directory = "plots", 
                           transformation.values = c(2,exp(1),3,4,5,6,7,8,9,10)){
  #logBases <- c(2,exp(1),3,4,5,6,7,8,9,10)
  #powers <- c(2,exp(1),3,4,5,6,7,8,9,10)
  #roots <- c(2,exp(1),3,4,5,6,7,8,9,10)
  
  record <- data.frame(
    index = index,
    feature = feature,
    values = data,
    flag = "Non-normal",
    transf = "",
    transf.value = 0
  )
  
  pvalues <- data.frame(matrix(vector(),3,length(transformation.values)))
  for (k in 1:ncol(pvalues)) {
    suppressWarnings({
      pvalues[1,k] <- shapiro.test(log(data,transformation.values[k]))[[2]]
      pvalues[2,k] <- shapiro.test(data^transformation.values[k])[[2]]
      pvalues[3,k] <- shapiro.test(data^(1/transformation.values[k]))[[2]]
    })
  }
  max.pvalue <- max(pvalues, na.rm = TRUE)
  max.pvalue.index <- which(pvalues == max.pvalue, arr.ind = TRUE)

  transf <- max.pvalue.index[1]
  transf.value.index <- max.pvalue.index[2]
  
  if(max.pvalue < 0.05) # Verify if a transformation normalized the data
    return(data.frame()) 
  
  if(transf == 1) { # Log transformation
    base <- transformation.values[transf.value.index]
    transformed <- log(data,base)
    if(base == exp(1))
      base <- "e"
    xlab <- paste0("$\\log_{",base,"}(",feature,")$")
    transformation <- paste0("LOG_",base)
    name.prefix <- paste0(plots.directory,"/HIST_",(index - offset),"_",transformation)
    generate.overlay.histogram(data,transformed,feature,name.prefix,xlab)
    record$values <- transformed
    record$transf <- "log"
    record$transf.value <- base
    return(record)
  }
  else if(transf == 2) { # Power transformation
    p <- transformation.values[transf.value.index]
    transformed <- data^p
    if(p == exp(1))
      p <- "e"
    xlab <- paste0("$(",feature,")^",p,"$")
    transformation <- paste0("POW_",p)
    name.prefix <- paste0(plots.directory,"/HIST_",(index - offset),"_",transformation)
    generate.overlay.histogram(data,transformed,feature,name.prefix,xlab)
    record$values <- transformed
    record$transf <- "power"
    record$transf.value <- p
    return(record)
  }
  else { # Root transformation
    r <- transformation.values[transf.value.index]
    transformed <- data^(1/r)
    if(r == exp(1))
      r <- "e"
    xlab <- paste0("$\\sqrt[",r,"]{",feature,"}$")
    transformation <- paste0("ROOT_",r)
    name.prefix <- paste0(plots.directory,"/HIST_",(index - offset),"_",transformation)
    generate.overlay.histogram(data,transformed,feature,name.prefix,xlab)
    record$values <- transformed
    record$transf <- "root"
    record$transf.value <- r
    return(record)
  }
}
