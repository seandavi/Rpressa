\name{readAgilent}
\Rdversion{1.1}
\alias{readAgilent}
\title{
  Read Agilent .txt files
}
\description{
  This is a fairly smart function for reading Agilent text files for
  miRNA, CGH, and one- and two-color expression arrays.
}
\usage{
readAgilent(fnames,
            path,
            dataColumns = NULL,
            annotationColumns = NULL,
            source = c("mirna", "cgh", "onecolorexpression",
                       "twocolorexpression"),
            exprs = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fnames}{
    A character() vector of file names
}
  \item{path}{
    The base file path for the associated file names.
}
  \item{dataColumns}{
    The column names of the data columns to be included.  A single
    matrix will be included in the assayData of the resulting
    \code{\linkS4class{AgilentSet}} object for each data column named
    here. Based on the source parameter, logical defaults are defined
    if this is NULL.
}
  \item{annotationColumns}{
    The column names of the annotation columns to be read.  These
    columns will be read from the first file and included in the
    featureData of the resulting \code{\linkS4class{AgilentSet}}
    object.  Based on the source parameter, logical defaults are defined
    if this is NULL.
}
  \item{source}{
    A single character string describing the type of array.  This
    selection defines the defaults for annotationColumns and dataColumns.
}
\item{exprs}{
  The integer index of the dataColumns vector to be assigned to the
  exprs() slot in the resulting \code{\linkS4class{AgilentSet}}.
}
  \item{\dots}{
    Passed to read.table.
}
}
\value{
  An object of class \code{\linkS4class{AgilentSet}}.
}
\author{
  Sean Davis <sdavis2@mail.nih.gov>
}
\seealso{
  \code{\linkS4class{AgilentSet}}
}
\examples{
}
\keyword{ IO }
