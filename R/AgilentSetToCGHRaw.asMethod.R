### Convert an AgilentSet to a CGHbase::cghRaw object
setAs("AgilentSet","cghRaw",function(from,to) {
  require(CGHbase)
  tmp = suppressMessages(splitAgilentChromLocs(fData(from)$SystematicName))
  colnames(tmp) = c('Chromosome','Start','End')
  tmp$Chromosome = suppressMessages(as.numeric(as.character(ChrNumeric(tmp$Chromosome))))
  ## Exclude unmapped probes and those that don't have
  ## the standard chromosome names (for human), such
  ## as chr6_random and chrM
  includedRows = !is.na(tmp$Chromosome)
  includedRows
  tmpfd = cbind(fData(from)[includedRows,],tmp[includedRows,])
  ret = new('cghRaw',phenoData=phenoData(from),
    experimentData=experimentData(from),
    annotation=annotation(from),
    copynumber=exprs(from)[includedRows,],
    featureData=as(tmpfd,"AnnotatedDataFrame"))
  ## Sort the cghRaw object before returning
  ret = ret[order(chromosomes(ret),bpstart(ret)),]
  return(ret)
})
