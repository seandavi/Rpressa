\name{feparams}
\Rdversion{1.1}
\alias{feparams}
\title{
  Access the Feature Extraction parameters from Agilent text file
}
\description{
    Access the Feature Extraction parameters from Agilent text file
}
\usage{
feparams(.Object)
}
\arguments{
  \item{.Object}{
    An object with an feparams slot, typically an "AgilentSet" class object.
  }
  }
\value{
  A data.frame with the feature extraction parameters in columns and
  samples in rows.
}
\author{
  Sean Davis <sdavis2@mail.nih.gov>
}
\seealso{
  \code{\link{festats}}, \code{\linkS4class{AgilentSet}}
}
\examples{

}
