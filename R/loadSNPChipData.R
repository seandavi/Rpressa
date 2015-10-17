loadSNPChipData <- function(fname,...) {
  f = gzfile(fname)
  dat <- read.table(f,sep="\t",header=TRUE,skip=9,nrow=100)
  colclasses <- sapply(dat,class)
  colclasses[-c(15:16,29:32)] <- "NULL"
  print(colclasses)
  f = gzfile(fname)
  dat2 <- read.table(f,sep="\t",header=TRUE,skip=9,colClasses=colclasses,...)
  return(dat2)
}
