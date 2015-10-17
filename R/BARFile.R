readParams <- function(con,n) {
  if(n>0) {
    for(i in 1:n) {
      parLength <- readBin(con,n=1,endian='big',what='integer')
      parName <- rawToChar(readBin(con,n=parLength,endian='big',what='raw'))
      parLength <- readBin(con,n=1,endian='big',what='integer')
      parVal <- rawToChar(readBin(con,n=parLength,endian='big',what='raw'))
    }
  }
}


readBar <- function(fname,verbose=TRUE) {
  f <- file(fname,'rb')
  x <- readBin(f,n=8,endian='big',what='raw')
  fileVersion <- readBin(f,n=1,endian='big',what='double')
  seek(f,12)
  numSeqs <- readBin(f,n=1,endian='big',what='integer')
  numDatCols <- readBin(f,n=1,endian='big',what='integer')
  fieldTypes <- readBin(f,n=numDatCols,endian='big',what='integer')
  numParPairs <- readBin(f,n=1,endian='big',what='integer')
  readParams(f,numParPairs)
  if(numSeqs>0) {
    print(numSeqs)
    for(i in 1:numSeqs) {
      x <- readBin(f,n=1,endian='big',what='integer')
      print(x)
      seqName <- rawToChar(readBin(f,n=x,endian='big',what='raw'))
      if(verbose) print(seqName)
      x <- readBin(f,n=1,endian='big',what='integer')
      print(x)
      seqGroup <- rawToChar(readBin(f,n=x,endian='big',what='raw'))
      if(verbose) print(seqGroup)
      x <- readBin(f,n=1,endian='big',what='integer')
      seqVersion <- rawToChar(readBin(f,n=x,endian='big',what='raw'))
      if(verbose) print(seqVersion)
      numParPairs <- readBin(f,n=1,endian='big',what='integer')
      readParams(f,numParPairs)
      x <- readBin(f,n=1,endian='big',what='integer')
      if(verbose) print(x)
      loc <- integer(x)
      val <- double(x)
      for(i in 1:x) {
        loc <- readBin(f,n=1,endian='big',what='integer')
        print(loc)
        val <- readBin(f,n=4,endian="big",what="raw")
        print(val)
      }
    }
  }
  close(f)
}
    
    
   
    
  
  
  
