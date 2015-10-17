\name{AgilentSet-class}
\Rdversion{1.1}
\docType{class}
\alias{AgilentSet-class}
\alias{initialize,AgilentSet-method}
\alias{as,AgilentSet-method}

\title{Class "AgilentSet"}
\description{ Derived from an eSet, adds two data.frames to the class
  for storing feparams and festats from the Agilent output}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AgilentSet", feparams, festats, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{feparams}:}{Object of class \code{"data.frame"}
      containing the output of the Agilent feature extraction
      parameters, one row per sample. }
    \item{\code{festats}:}{Object of class \code{"data.frame"}
      containing the output of the Agilent feature extraction stats, one
    row per sample.}
    \item{\code{assayData}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{ExpressionSet}"}, directly.
Class \code{"\linkS4class{eSet}"}, by class "ExpressionSet", distance 2.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "ExpressionSet", distance 3.
Class \code{"\linkS4class{Versioned}"}, by class "ExpressionSet", distance 4.
}
\section{Methods}{
  \describe{
    \item{feparams}{\code{signature(.Object = "AgilentSet")}: accessor
      for the feparams slot }
    \item{festats}{\code{signature(.Object = "AgilentSet")}: accessor
      for the festats slot }
    \item{initialize}{\code{signature(.Object = "AgilentSet")}: ... }
	 }
    \item{as}{\code{signature(from=AgilentSet,to="cghRaw")}: Coerce an AgilentSet to 
      CGHbase::cghRaw.}
	 }
}
\author{ Sean Davis <sdavis2@mail.nih.gov> }
\seealso{
  \code{\link{readAgilent}} and \code{\linkS4class{eSet}} 
}
\examples{
showClass("AgilentSet")
}
\keyword{classes}
