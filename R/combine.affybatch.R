setMethod("combine",signature=c("AffyBatch","AffyBatch"),function(x,y,...) {
    warning("this is not a robust function right now.  Data could come in in different orders and this function would not catch it")
      if(annotation(x)!=annotation(y) |
              cdfName(x)!=cdfName(y) |
              nrow(x)!=nrow(y) |
              ncol(x)!=ncol(y)) {
            stop("combine works only with datasets from the same chip type")
          }
      experimentdata <- combine(experimentData(x),experimentData(y))
      phenodata <- combine(phenoData(x),phenoData(y))
      featuredata <- combine(featureData(x),featureData(y))
      protocoldata <- combine(protocolData(x),protocolData(y))
      cat('combining exprs...\n')
      f1 <- match(featureNames(x),featureNames(y))
      if(any(is.na(f1))) {
            stop("Featurenames did not match between affybatches!")
          }
      exps <- cbind(exprs(x),exprs(y))
      cat('creating affybatch...\n')
      tmp <- new("AffyBatch",cdfName=cdfName(x),nrow=nrow(x),
                              ncol=ncol(x),annotation=annotation(x),
                              phenoData=phenodata,featureData=featuredata,
                              protocolData=protocoldata,experimentData=experimentdata,
                              assayData=assayDataNew(exprs=exps))
      return(tmp)
  })
