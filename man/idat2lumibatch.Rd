\name{idat2lumibatch}
\alias{idat2lumibatch}
\title{Read idat files into a lumiBatch object}
\usage{
idat2lumibatch(filenames, annotation)
}
\arguments{
  \item{filenames}{A character vector of the full path for
  each idat file to be read}

  \item{annotation}{The name, as a character string (eg.,
  "illuminaHumanv4.db"), of the annotation package to be
  used for identifying bead types and for inserting feature
  data.}
}
\description{
Given a set of idat filenames and an annotation package,
returns a lumiBatch object with detection p-values
calculated.
}

