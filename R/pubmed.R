pubmedQuery = function(search) {
  library(XML)
   return(xmlRoot(xmlTreeParse(sprintf("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=%s&tool=R",search))))
}

getQueryCounts = function(genes,terms,baseline=NULL) {
  retmat = matrix(rep(0,length(terms)*length(genes)),nc=length(terms))
  colnames(retmat)=terms
  rownames(retmat)=genes
  pb = txtProgressBar(max=length(terms)*length(genes),style=3)

  qcounts = function(doc) {
      if(length(getNodeSet(doc,'//ErrorList'))>0) {
                                        # If we are here, an error or "not found" occurred.
          return(0)
      } else {
          return(as.integer(xmlValue(getNodeSet(doc,"/eSearchResult/Count")[[1]])))
      }
  }
      

  i = 0
  for(gene in genes) {
    for(term in terms) {
      i=i+1
      setTxtProgressBar(pb,i)
      search = sprintf('"%s" AND "%s"',term,gene)
      if(!is.null(baseline)) {
        search = sprintf('%s AND "%s"',search,baseline)
      }
      doc=pubmedQuery(search)
      retmat[gene,term]=qcounts(doc)
    }
    Sys.sleep(0.4)
  }
  genecounts = sapply(genes,function(x) {
      search = sprintf('"%s"',x)
      if(!is.null(baseline)) {
        search = sprintf('%s AND "%s"',search,baseline)
      }
      Sys.sleep(0.4)
      return(qcounts(pubmedQuery(search)))
  })
  names(genecounts)=genes
  termcounts = sapply(terms,function(x) {
      search = sprintf('"%s"',x)
      if(!is.null(baseline)) {
        search = sprintf('%s AND "%s"',search,baseline)
      }
      Sys.sleep(0.4)
      return(qcounts(pubmedQuery(search)))
  })
  names(termcounts)=terms
  return(list(genecounts=genecounts,
              termcounts=termcounts,
              matrix=retmat))
}

checkCounts = function(querycountoutput) {
    x = querycountoutput
    res = sapply(names(x$genecounts),function(gene) {
        mapply(x$matrix[gene,],x$termcounts,x$genecounts[gene],FUN=function(a,b,c) return(fisher.test(matrix(c(a,c,b,21000000),nc=2))$p.value))})
    colnames(res) = names(x$genecounts)
    return(t(res))
}
