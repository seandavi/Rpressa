require(Biobase)
setClass('AgilentSet',
         representation(feparams="data.frame",  ### Need to convert over to annotatedDataFrame
                        festats="data.frame"),  ### Need to convert over to annotatedDataFrame
         contains="ExpressionSet")

setMethod('initialize','AgilentSet',function(.Object,feparams,festats,...) {
  .Object@feparams=feparams
  .Object@festats=festats
  callNextMethod(.Object,...)
})

setGeneric('feparams',function(.Object) {standardGeneric('feparams')})
setMethod('feparams','AgilentSet',function(.Object) {.Object@feparams})

setGeneric('festats',function(.Object) {standardGeneric('festats')})
setMethod('festats','AgilentSet',function(.Object) {.Object@festats})

convertDataType <- function(AgilentType,val) {
  convertfcn <- switch(AgilentType,
                     integer=as.integer,
                     float=as.numeric,
                     boolean=as.logical,
                     text=as.character)
  return(convertfcn(val))
}
                     


readAgilentHeader <- function(con) {
  dat <- readLines(con,n=7)
  dattypes <- as.list(strsplit(dat[1],"\t")[[1]][-1])
  feparams <- as.list(strsplit(dat[3],"\t")[[1]][-1])
  feparams <- mapply(convertDataType,dattypes,feparams,SIMPLIFY=FALSE)
  names(feparams) <- unlist(strsplit(dat[2],"\t")[[1]][-1])
  dattypes <- as.list(strsplit(dat[5],"\t")[[1]][-1])
  festats <- as.list(strsplit(dat[7],"\t")[[1]][-1])
  festats <- mapply(convertDataType,dattypes,festats,SIMPLIFY=FALSE)
  names(festats) <- unlist(strsplit(dat[6],"\t")[[1]][-1])
  return(list(feparams=feparams,festats=festats))
}

fastReadTable <- function(con,columns,ntest=100,...) {
  loc <- seek(con)
  dat <- read.table(con,nrows=ntest,...)
  colclasses <- rep(NA,ncol(dat))
  cols <- match(columns,colnames(dat))
  if(any(is.na(cols))) {
    warning(paste("Column(s)",paste(columns[is.na(cols)],sep=", "),"missing in data file"))
  }
  cols <- cols[!is.na(cols)]
  colclasses[-cols] <- "NULL"
  seek(con,loc)
  newdat <- read.table(con,colClasses=colclasses,...)
}
  
  

readAgilent <- function(fnames,path,dataColumns=NULL,annotationColumns=NULL,source=c('mirna','cgh','onecolorexpression','twocolorexpression'),exprs=1,...) {
  if(is.null(dataColumns)) {
    dataColumns <- switch(source,
                          mirna=c('gProcessedSignal','gIsPosAndSignif','gTotalProbeSignal','gTotalProbeError','gIsGeneDetected','gBGMedianSignal'),
                          cgh=c('LogRatio','LogRatioError','rProcessedSignal','gProcessedSignal'),
                          onecolorexpression=c('gProcessedSignal','gProcessedSignalError','gIsPositiveAndSignif','gProcessedBackground'),
                          twocolorexpression=c('gProcessedSignal','gProcessedSignalError','gIsPositiveAndSignif','gProcessedBackground','rProcessedSignal','rProcessedSignalError','rIsPositiveAndSignif','rProcessedBackground'))
  }
  if(is.null(annotationColumns)) {
    annotationColumns <- switch(source,
                                mirna=c('Row','Col','chr_coord','ProbeUID','ControlType','ProbeName','GeneName','SystematicName','Description'),
                                twocolorexpression=c('Row','Col','ProbeUID','ControlType','ProbeName','GeneName','SystematicName','Description'),
                                cgh=c('Row','Col','ProbeUID','ControlType','ProbeName','SystematicName'),
                                onecolorexpression=c('Row','Col','ProbeUID','ControlType','ProbeName','GeneName','SystematicName','Description'))
  }
  firsttime <- TRUE
  annotdat <- NULL
  signaldat <- list()
  feparamslist <- list()
  festatslist <- list()
  for(i in 1:length(fnames)) {
    print(fnames[i])
    fullname <- fname <- fnames[i]
    if(!is.null(path)) {
      fullname <- file.path(path,fname)
    }
    con <- file(fullname,'r')
    columns <- ""
    if(firsttime) {
      columns <- union(dataColumns,annotationColumns)
    } else {
      columns <- dataColumns
    }
    headerinfo <- readAgilentHeader(con)
    feparamslist[[i]] <- headerinfo$feparams
    festatslist[[i]] <- headerinfo$festats
    
    dat <- fastReadTable(con,columns,ntest=1000,skip=2,sep="\t",header=TRUE,quote="",comment.char="",...)
    if(firsttime) {
      tmpdatcols <- match(annotationColumns,colnames(dat))
      tmpdatcols2 <- tmpdatcols[!is.na(tmpdatcols)]
      annotdat <- dat[,tmpdatcols2]
      colnames(annotdat) <- annotationColumns[!is.na(tmpdatcols)]
#      annotdat <- switch(source,
#                         mirna=cbind(annotdat,splitAgilentChromLocs(dat$chr_coord)),
#                         cgh=cbind(annotdat,splitAgilentChromLocs(dat$SystematicName)),
#                         annotdat)
      tmpdatcols <- match(dataColumns,colnames(dat))
      tmpdatcols <- tmpdatcols[!is.na(tmpdatcols)]
      dat <- dat[,tmpdatcols]
      firsttime <- FALSE
    }
    signaldat[[i]] <- dat
    if(length(unique(sapply(signaldat,nrow)))!=1) {
      stop("Most recent file did not have same number of rows as others")
    }
    close(con)
  }
  a <- new.env(hash=FALSE)
  for(i in 1:length(dataColumns)) {
    tmp <- as.matrix(sapply(signaldat,function(x) {x[,i]}))
    colnames(tmp) <- fnames
    if(i==exprs) {
      assign('exprs',tmp,envir=a)
    } else {
      assign(dataColumns[i],tmp,envir=a)
    }
  }
  if(is.null(annotdat)) {
    return(new("AgilentSet",assayData=a,
               feparams=data.frame(do.call(rbind,feparamslist)),
               festats=data.frame(do.call(rbind,festatslist))))
  } else {
    if("SystematicName" %in% colnames(annotdat) & source=='cgh') {
      annotdat <- cbind(annotdat,splitAgilentChromLocs(annotdat$SystematicName))
    }
    return(new("AgilentSet",assayData=a,featureData=as(annotdat,"AnnotatedDataFrame"),
               feparams=data.frame(do.call(rbind,feparamslist)),
               festats=data.frame(do.call(rbind,festatslist))))
  }
}

splitAgilentChromLocs <- function(systematicName) {
  tmp <- gsub('[:-]',':',as.character(systematicName))
  tmp2 <- data.frame(do.call(rbind,strsplit(as.character(tmp),':')))
  colnames(tmp2) <- c('chrom','start','end')
  tmp2[,2]=as.integer(as.character(tmp2[,2]))
  tmp2[,3]=as.integer(as.character(tmp2[,3]))
  return(tmp2)
}
