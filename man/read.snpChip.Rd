\name{read.snpChip}
\Rdversion{1.1}
\alias{read.snpChip}
\title{
  Read in split Illumina snp chip files
}
\description{
  Reads the SNP_Map.txt, Sample_Map.txt, and the one-sample-per-file
  export files from the GenomeStudio export.
}
\usage{
read.snpChip(path='.',file.prefix,sampleMapFile=file.path(path,"Sample_Map.txt"),
           snpMapFile=file.path(path,"SNP_Map.txt"),...) {
}
}
\arguments{
  \item{path}{
    The string representing the path to the GenomeStudio export files.
}
  \item{file.prefix}{
    The GenomeStudio export file prefix that will be prepended to the
file names for each sample.  In the end, the file names will be of the
form "file.prefix##.txt" where "##" is the sample number.
}
  \item{sampleMapFile}{
    The default will specify the default file name from the GenomeStudio
output in the same directory as the exported data files.  This file
contains the sample annotations.  The resulting data frame becomes the
phenoData for the ExpressionSet.
}
  \item{SNPMapFile}{
    The default will specify the default file name from the GenomeStudio
output in the same directory as the exported data files.  This file
contains the SNP annotations and will end up in the featureData for the
resulting ExpressionSet.
}
  \item{\dots}{
    Passed to read.table()
}
}
\details{
  
}
\value{
An "ExpressionSet" object.
}
\author{
  Sean Davis <seandavi@gmail.com>
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
\keyword{ IO }
