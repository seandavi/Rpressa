oligonucleotideFrequencyByChromosome <- function(chromosome,
                                                 start,
                                                 end,
                                                 BSgenomeObject,
                                                 width=1) {
  require("Biostrings")
  require("BSgenome")
  uniqChroms <- unique(chromosome)
  firsttime <- TRUE
  dat <- NULL
  start[start<1] <- 1
  for(i in uniqChroms) {
    cat(i,'\n')
    if(!(i %in% names(BSgenomeObject))) {
      next
    }
    idx <- which(chromosome==i)
    dnastring <- DNAString(BSgenomeObject[[i]])
    end[idx][end[idx]>length(dnastring)] <- length(dnastring)
    v <- Views(dnastring,start=start[idx],end=end[idx])
    tmp <- oligonucleotideFrequency(v,width)
    if(firsttime) {
      dat <- matrix(nrow=length(chromosome),ncol=ncol(tmp))
      colnames(dat) <- colnames(tmp)
      firsttime <- FALSE
    }
    dat[idx,] <- tmp
  }
  return(dat)
}
