#' Read idat files into a lumiBatch object
#'
#' Given a set of idat filenames and an annotation package,
#' returns a lumiBatch object with detection p-values calculated.
#'
#' @param filenames A character vector of the full path for each idat file to be read
#' @param annotation The name, as a character string (eg., "illuminaHumanv4.db"),
#' of the annotation package to be used for identifying bead types and for
#' inserting feature data.
#'
#' @export
#' 
idat2lumibatch <- function(filenames,annotation) {
  require(illuminaio)
  require(lumi)
  require(annotation,character.only=TRUE)
  idatlist = lapply(filenames,readIDAT)
  exprs = sapply(idatlist,function(x) {
    return(x$Quants$MeanBinData)})
  se.exprs = sapply(idatlist,function(x) {
    return(x$Quants$DevBinData/sqrt(x$Quants$NumGoodBeadsBinData))})
  beadNum = sapply(idatlist,function(x) {
    return(x$Quants$NumGoodBeadsBinData)})
  rownames(exprs)=rownames(se.exprs)=rownames(beadNum)=idatlist[[1]]$Quants$CodesBinData
  colnames(exprs)=colnames(se.exprs)=colnames(beadNum)=sapply(idatlist,function(x) {
    return(paste(x$Barcode,x$Section,sep="_"))})
  pd = data.frame(Sentrix=colnames(exprs))
  rownames(pd)=colnames(exprs)
  conn = illuminaHumanv4_dbconn()
  tmp = dbGetQuery(conn,'select * from ExtraInfo')
  lb = new("LumiBatch",exprs=exprs,se.exprs=se.exprs,beadNum=beadNum,
    phenoData=AnnotatedDataFrame(pd))
  z = match(featureNames(lb),tmp$ArrayAddress)
  lb = lb[!is.na(z),]
  z = match(featureNames(lb),tmp$ArrayAddress)
  fData(lb)=tmp[z,]
  negprobes = which(fData(lb)$ReporterGroupName=='negative')
  browser()
  assayDataElement(lb,'pvals') = apply(exprs(lb),2,function(a) {
      return(illuminaPvalCalculation(a,negprobes))})
  return(lb)
}


#' Encapsulates a vectorized p-value calculation for Illumina microarrays
#'
#' The Illumina expression arrays use negative control probes for
#' determining the p-value of expression.  Refer to the Genome Studio
#' Gene Expression Module documentation for details.
#' 
#' @param values A vector of numeric values associated with a single sample
#' @param negProbeIdx An index defining the values that represent the negative
#' control probes
#'
#' @export
#' 
illuminaPvalCalculation = function(values,negProbeIdx) {
    ec = ecdf(values[negProbeIdx])
    return(1-ec(values))
}
