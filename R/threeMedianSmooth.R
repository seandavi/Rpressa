threeMedianSmooth <- function(x) {
  tmp <- rowMedians(embed(x,3))
  return(c(mean(c(x[1],x[2])),tmp,mean(c(x[length(x)-1],x[length(x)]))))
}
  
