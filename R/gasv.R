generateRangesForNormal = function(dat) {
  #
  # used for output file from GASV
  #   cluster format:
  # 
  #
  #   X.Cluster_ID. LeftChr. LeftBreakPoint. RightChr. RightBreakPoint. Num.PRS. Localization. Type.
  #   1           c66        1   964969,965167         1    965255,965453        5         123.2     D
  #   2           c70        1   965872,966311         1    966242,966646        7          -1.0     V
  #   3          c217        1 1531128,1531752         1  1531567,1532173        4          -1.0     D
  #   4          c254        1 1649585,1649931         1  1650277,1650467        8          -1.0     D
  #   5          c490        1 2585206,2585510         1  2585287,2585591       10         210.0     D
  #   6          c495        1 2586111,2586399         1  2629238,2629526        4         194.8     V
  
  pos = do.call('rbind',strsplit(c(as.character(dat[,3]),as.character(dat[,5])),','))
  seqnames = c(as.character(dat[,2]),as.character(dat[,4]))
  seqnames[seqnames=='23']='X'
  seqnames[seqnames=='24']='Y'
  return(GRanges(seqnames=seqnames,ranges=IRanges(start=as.numeric(as.character(pos[,1])),
                                                  end=as.numeric(as.character(pos[,2])))))
}

generateRangeForTumor = function(dat) {
  #
  # used for output file from GASV
  #   cluster format:
  # 
  #
  #   X.Cluster_ID. LeftChr. LeftBreakPoint. RightChr. RightBreakPoint. Num.PRS. Localization. Type.
  #   1           c66        1   964969,965167         1    965255,965453        5         123.2     D
  #   2           c70        1   965872,966311         1    966242,966646        7          -1.0     V
  #   3          c217        1 1531128,1531752         1  1531567,1532173        4          -1.0     D
  #   4          c254        1 1649585,1649931         1  1650277,1650467        8          -1.0     D
  #   5          c490        1 2585206,2585510         1  2585287,2585591       10         210.0     D
  #   6          c495        1 2586111,2586399         1  2629238,2629526        4         194.8     V
  posl = do.call('rbind',strsplit(as.character(dat[,3]),','))
  posr = do.call('rbind',strsplit(as.character(dat[,5]),','))
  seqnamesl = as.character(dat[,2])
  seqnamesl[seqnamesl=='23']='X'
  seqnamesl[seqnamesl=='24']='Y'
  grl=GRanges(seqnames=seqnamesl,ranges=IRanges(start=as.numeric(as.character(posl[,1])),
                                                end=as.numeric(as.character(posl[,2]))))
  seqnamesr = as.character(dat[,4])
  seqnamesr[seqnamesr=='23']='X'
  seqnamesr[seqnamesr=='24']='Y'
  grr=GRanges(seqnames=seqnamesr,ranges=IRanges(start=as.numeric(as.character(posr[,1])),
                                                end=as.numeric(as.character(posr[,2]))))
  return(list(left=grl,right=grr))
}

