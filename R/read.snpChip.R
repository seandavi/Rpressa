readSNPMap <-
  function(fname,...) {
    dat <- read.table(fname,sep="\t",header=TRUE,...)
    return(dat)
  }

readSampleMap <-
  function(fname,...) {
    dat <- read.table(fname,sep="\t",header=TRUE)
    rownames(dat) <- dat$Index
    return(dat)
  }
    

read.snpChip <-
  function(path='.',file.prefix,sampleMapFile=file.path(path,"Sample_Map.txt"),
           snpMapFile=file.path(path,"SNP_Map.txt"),...) {
    require(Biobase)
    colclasses <- rep("NULL",34)
    colclasses[29:32] <- "numeric"
    colclasses[15:16] <- "factor"
    samps <- readSampleMap(sampleMapFile)
    annot <- readSNPMap(snpMapFile)
    a <- new("ExpressionSet")
    featureData(a) <- as(annot,"AnnotatedDataFrame")
    phenoData(a) <- as(samps,"AnnotatedDataFrame")
    x <- file.path(path,paste(file.prefix,1:nrow(samps),'.txt',sep=""))
    d <- list()
    tmp <- lapply(x,function(y) {
      print(y)
      tmp1 <- read.table(y,colClasses=colclasses,sep="\t",header=TRUE,skip=10,...)
      return(tmp1)
    })
    d <- list()
    for(i in 1:ncol(tmp[[1]])) {
      d[[i]] <- do.call(cbind,lapply(tmp,'[',i))
    }
    names(d) <- c('Allele1','Allele2','XRaw','YRaw','baf','exprs')
    assayData(a) <- assayDataNew(exprs=d$exprs,baf=d$baf,xraw=d$XRaw,
                                 yraw=d$YRaw,allele1=d$Allele1,
                                 allele2=d$Allele2)
    sampleNames(a) <- 1:length(x)
    featureNames(a) <- annot[,1]
    return(a)
  }
    

