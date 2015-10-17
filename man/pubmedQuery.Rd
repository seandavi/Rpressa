\name{pubmedQuery}
\alias{pubmedQuery}
\title{
Return results of a PubMed query
}
\description{
Given a search string, perform a PubMed query and return the XML result.
}
\usage{
pubmedQuery(search)
}
\arguments{
  \item{search}{
  A single seach string as entered in PubMed.
  }
}
\value{
The raw XML result from an eSearch call to PubMed.
}
\references{
See \url{http://www.ncbi.nlm.nih.gov/books/NBK25500/} for details.
}
\author{
Sean Davis <sdavis2@mail.nih.gov>
}
\seealso{
\code{\link{getQueryCounts}}
}
\examples{
search = "science[journal] AND breast cancer AND 2008[pdat]"
res = pubmedQuery(search)
res
}
\keyword{ misc }
