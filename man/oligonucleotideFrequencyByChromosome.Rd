\name{oligonucleotideFrequencyByChromosome}
\Rdversion{1.1}
\alias{oligonucleotideFrequencyByChromosome}
\title{
  Find oligonucleotide frequencies for multiple chromosomes
}
\description{
  The oligonucleotideFrequency function is useful for finding things
  like GC percentage, etc., from a DNAString object (or others).
  However, if one wants to use the combination of BSgenome and
  oligonucleotideFrequency to find these parameters over the genome,
  this is a tedious process.  This little function streamlines the
  process for doing things like finding the GC percentage of probes or
  regions across the entire genome of interest.
}
\usage{
oligonucleotideFrequencyByChromosome(chromosome, start, end, BSgenomeObject, width = 1)
}
\arguments{
  \item{chromosome}{
    A character vector of the chromosome of each region.  The vector
    members need to match the names used in the BSgenome data package
    (typically, these are going to look like the UCSC chromosomes,
    "chr1", "chr2", ....)
}
  \item{start}{
    An integer vector (1-based) of the start locations for the
    oligonucleotide frequency calculations.  Start locations that, for
    whatever reason, are less than 1 are set to 1.
}
  \item{end}{
    An integer vector (1-based) of the end locations for the
    oligonucleotide frequency calculations.  End locations that, for
    whatever reason, are greater than the length of the respective
    chromosome are set to the length of the chromosome.
}
  \item{BSgenomeObject}{
    The BSgenome object.  For example, if the genome of interest is
    drawn from the "BSgenome.Hsapiens.UCSC.hg18", then the variable that
    holds the genome sequences is Hsapiens.  It is this actual variable
    that should be used here.
}
  \item{width}{
    Integer value describing the width of the nmers to check (1 for A,
    C, T, G; 2 for dinucleotides; etc.)
}
}
\value{
  An integer matrix with the counts of oligonucleotides, one row for
  each location, in the same order as the input chromosome, start, end.
}
\author{
  Sean Davis <sdavis2@mail.nih.gov>
}
\seealso{
  \code{\link{oligonucleotideFrequency}}
}
\examples{
require(BSgenome)
require(BSgenome.Hsapiens.UCSC.hg18)
oligonucleotideFrequencyByChromosome(chromosome=c('chr1','chr2','chr3'),
                                     start=c(1,1,1),
                                     end=c(10000,10000,9000),
                                     BSgenomeObject=Hsapiens,
                                     width=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
