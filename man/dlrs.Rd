\name{dlrs}
\Rdversion{1.1}
\alias{dlrs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate the Derivative Log Root Square (DLRS)
}
\description{
  The DLRS is a robust estimate of the standard deviation.  It is useful
  for determining the local spread of a distribution, such as for
  determining the noise in a CGH array.  
}
\usage{
dlrs(x)
}
\arguments{
  \item{x}{A numeric vector.
}
}
\details{
  The DLRS is a robust estimate of the standard deviation.  It is useful
  for determining the local spread of a distribution, such as for
  determining the noise in a CGH array.    
}
\value{
  A single numeric value representing the DLRS.
}
\author{
  Sean Davis <seandavi@gmail.com>
}
\seealso{
  \code{\link{mad}}
}
\examples{
x <- rnorm(10000)
dlrs(x)
sd(x)
# note that a shift of a large chunk of the data
# does not significantly alter the dlrs, but it
# does alter the sd.
x[5000:7000] <- x[5000:7000]+10
dlrs(x)
sd(x)
}
\keyword{ robust }
\keyword{ distribution }
