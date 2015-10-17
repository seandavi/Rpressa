\name{targeted}
\alias{targeted}
\docType{data}
\title{
Sample AlignedRead data set for targeted sequencing application
}
\description{
This is a dataset that was generated using in-solution capture technology to capture some exons of genes of interest.  The regions of capture are described in another related dataset.
}
\usage{data(targeted)}
\format{
  The format is:
Formal class 'AlignedRead' [package "ShortRead"] with 8 slots
  ..@ chromosome  : Factor w/ 66767 levels "0:0:1","0:0:10",..: 66767 66767 22181 7859 66767 66767 66767 66767 66767 66767 ...
  ..@ position    : int [1:6349344] NA NA NA NA NA NA NA NA NA NA ...
  ..@ strand      : Factor w/ 3 levels "-","+","*": NA NA NA NA NA NA NA NA NA NA ...
  ..@ alignQuality:Formal class 'NumericQuality' [package "ShortRead"] with 1 slots
  .. .. ..@ quality: int [1:6349344] 0 0 0 0 0 0 0 0 0 0 ...
  ..@ alignData   :Formal class 'AlignedDataFrame' [package "ShortRead"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	7 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr [1:7] "Analysis pipeline run" "Flow cell lane" "Flow cell tile" "Cluster x-coordinate" ...
  .. .. ..@ data             :'data.frame':	6349344 obs. of  7 variables:
  .. .. .. ..$ run      : Factor w/ 1 level "90814": 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. ..$ lane     : int [1:6349344] 5 5 5 5 5 5 5 5 5 5 ...
  .. .. .. ..$ tile     : int [1:6349344] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. ..$ x        : int [1:6349344] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. .. ..$ y        : int [1:6349344] 1947 1976 623 330 552 90 133 315 410 447 ...
  .. .. .. ..$ filtering: Factor w/ 2 levels "Y","N": 2 2 2 2 2 2 2 2 2 2 ...
  .. .. .. ..$ contig   : Factor w/ 1 level "": 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..@ dimLabels        : chr [1:2] "readName" "alignColumn"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ quality     :Formal class 'SFastqQuality' [package "ShortRead"] with 1 slots
  .. .. ..@ quality:Formal class 'BStringSet' [package "Biostrings"] with 5 slots
  .. .. .. .. ..@ super          :Formal class 'BString' [package "Biostrings"] with 6 slots
  .. .. .. .. .. .. ..@ xdata          :Formal class 'RawPtr' [package "IRanges"] with 2 slots
  .. .. .. .. .. .. .. .. ..@ xp                    :<externalptr> 
  .. .. .. .. .. .. .. .. ..@ .link_to_cached_object:<environment: 0x566412c> 
  .. .. .. .. .. .. ..@ offset         : int 0
  .. .. .. .. .. .. ..@ length         : int 253973760
  .. .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. .. ..@ elementType    : chr "ANYTHING"
  .. .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots
  .. .. .. .. .. .. ..@ start          : int [1:6349344] 1 41 81 121 161 201 241 281 321 361 ...
  .. .. .. .. .. .. ..@ width          : int [1:6349344] 40 40 40 40 40 40 40 40 40 40 ...
  .. .. .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. .. .. ..@ metadata       : list()
  .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. ..@ elementType    : chr "ANYTHING"
  .. .. .. .. ..@ metadata       : list()
  ..@ sread       :Formal class 'DNAStringSet' [package "Biostrings"] with 5 slots
  .. .. ..@ super          :Formal class 'DNAString' [package "Biostrings"] with 6 slots
  .. .. .. .. ..@ xdata          :Formal class 'RawPtr' [package "IRanges"] with 2 slots
  .. .. .. .. .. .. ..@ xp                    :<externalptr> 
  .. .. .. .. .. .. ..@ .link_to_cached_object:<environment: 0x566412c> 
  .. .. .. .. ..@ offset         : int 0
  .. .. .. .. ..@ length         : int 253973760
  .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. ..@ elementType    : chr "ANYTHING"
  .. .. .. .. ..@ metadata       : list()
  .. .. ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots
  .. .. .. .. ..@ start          : int [1:6349344] 1 41 81 121 161 201 241 281 321 361 ...
  .. .. .. .. ..@ width          : int [1:6349344] 40 40 40 40 40 40 40 40 40 40 ...
  .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. ..@ metadata       : list()
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ elementType    : chr "ANYTHING"
  .. .. ..@ metadata       : list()
  ..@ id          :Formal class 'BStringSet' [package "Biostrings"] with 5 slots
  .. .. ..@ super          :Formal class 'BString' [package "Biostrings"] with 6 slots
  .. .. .. .. ..@ xdata          :Formal class 'RawPtr' [package "IRanges"] with 2 slots
  .. .. .. .. .. .. ..@ xp                    :<externalptr> 
  .. .. .. .. .. .. ..@ .link_to_cached_object:<environment: 0x566412c> 
  .. .. .. .. ..@ offset         : int 0
  .. .. .. .. ..@ length         : int 0
  .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. ..@ elementType    : chr "ANYTHING"
  .. .. .. .. ..@ metadata       : list()
  .. .. ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots
  .. .. .. .. ..@ start          : int [1:6349344] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. .. ..@ width          : int [1:6349344] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. .. .. ..@ NAMES          : NULL
  .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. ..@ elementType    : chr "integer"
  .. .. .. .. ..@ metadata       : list()
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ elementType    : chr "ANYTHING"
  .. .. ..@ metadata       : list()
}
\source{
The data were generated in the Genetics Branch, National Cancer Institute, National Institutes of Health by Ogan Abaan.
}

\examples{
data(targeted)
targeted
}
\keyword{datasets}
