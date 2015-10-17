setGeneric('viewMeans',function(x,na.rm=FALSE) standardGeneric('viewMeans'))
setMethod('viewMeans','RleViews',function(x,na.rm=FALSE) {
  return(viewSums(x,na.rm=na.rm)/width(x))
})
setMethod('viewMeans','RleViewsList',function(x,na.rm=FALSE) {
  return(NumericList(sapply(x,viewMeans,na.rm=na.rm)))
})
setMethod('viewMeans','XIntegerViews',function(x,na.rm=FALSE) {
  return(viewSums(x,na.rm=na.rm)/width(x))
})
