\name{threeMedianSmooth}
\Rdversion{1.1}
\alias{threeMedianSmooth}
\title{
  Moving three-point median smoothing
}
\description{
  A three-point moving median is a very effective way of robustly
  removing outlier single points without removing effects in a vector
  that persist over more than one point.
}
\usage{
threeMedianSmooth(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector sorted in the order of interest (for
    example, chromosome order)}
}
\details{
  This is pretty simple.  The data are embedded with length=3 and then
  rowMedians is called.
}
\value{
  A numeric vector the same length as the input vector.
}
\author{
  Sean Davis <sdavis2@mail.nih.gov>
}
\examples{
x <- rnorm(100)
x[4] <- 10
summary(x)
y <- threeMedianSmooth(x)
summary(y)
}
\keyword{ stats }
