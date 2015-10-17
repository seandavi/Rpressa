readSailfish = function(filenames = dir('.',pattern='quant.sf',recursive=TRUE),
    samplenames=sub('/quant.sf','',filenames)) {
    d = sapply(filenames,read.delim,skip=4,simplify=FALSE)
    tpm = do.call(cbind,lapply(d,function(x) {x[,3]}))
    rpkm = do.call(cbind,lapply(d,function(x) {x[,4]}))
    colnames(tpm)=samplenames
    colnames(rpkm)=samplenames
    genes = d[[1]][,1:2]
    colnames(genes)=c('Transcript','Length')
    ev = new.env()
    assign('exprs',tpm,ev)
    assign('rpkm',rpkm,ev)
    eset = ExpressionSet(assayData=ev,featureData=as(genes,'AnnotatedDataFrame'))
    return(eset)
}
