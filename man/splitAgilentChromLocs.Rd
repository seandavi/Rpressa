\name{splitAgilentChromLocs}
\Rdversion{1.1}
\alias{splitAgilentChromLocs}
\title{
  Split Chromosome Locations from Agilent CGH arrays
}
\description{
  Agilent records the chromosome location information into a single
  column, typically called SystematicName.  This function just splits
  that column into chromosome, start, and end.
}
\usage{
splitAgilentChromLocs(systematicName)
}
\arguments{
  \item{systematicName}{
    A column from the Agilent array cgh data, typically with name
    SystematicName that contains the chromosome location information.
}
}
\value{
  A three-column data.frame with names "chrom", "start", and "end".
  Note that there will be NA values when the probe represents a control
  probe on the array.
}
\author{
  Sean Davis <sdavis2@mail.nih.gov>
}
\examples{
}
