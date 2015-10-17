\name{illuminaPvalCalculation}
\alias{illuminaPvalCalculation}
\title{Encapsulates a vectorized p-value calculation for Illumina microarrays}
\usage{
illuminaPvalCalculation(values, negProbeIdx)
}
\arguments{
  \item{values}{A vector of numeric values associated with
  a single sample}

  \item{negProbeIdx}{An index defining the values that
  represent the negative control probes}
}
\description{
The Illumina expression arrays use negative control probes
for determining the p-value of expression.  Refer to the
Genome Studio Gene Expression Module documentation for
details.
}

