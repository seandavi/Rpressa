\name{getQueryCounts}
\alias{getQueryCounts}
\title{
Return query counts for batch PubMed queries
}
\description{
Trying to discover the biological meaning of a set of genes can
be challenging.  This function simply does batch PubMed queries
and returns a matrix of record counts in matrix form with genes
as rows and terms as columns.
}
\usage{
getQueryCounts(genes, terms, baseline=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genes}{
    A character vector of gene names (or miRNAs, etc.).  These will
    become the rows of the query count matrix returned.
  }
  \item{terms}{
    A character vector of terms to pair with each gene.  These will
    become the columns of the query count matrix returned.
  }
  \item{baseline}{
    A string to include with ALL queries.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A list with three items:
  
A vector of the number of of pubmed results for each gene, including the baseline terms if specified

A vector of the number of of pubmed results for each term, including the
baseline terms if specified

A numeric matrix with column names taken from "terms" and 
row names taken from "genes".  Each entry in the matrix
is the count of records in PubMed for that pair of gene/term.
}
\references{
\url{http://www.ncbi.nlm.nih.gov/pubmed/}
}
\author{
  Sean Davis <sdavis2@mail.nih.gov>
}
\seealso{
\code{\link{pubmedQuery}}
}
\examples{
genes = c('BRCA1','AKT','PIK3CA','CDKN2A')
terms = c('cancer','cell cycle','tumor suppressor','oncogene')
qcounts = getQueryCounts(genes,terms)
qcounts
}
\keyword{ misc }
