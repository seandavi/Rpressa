mantaVCF2GRanges = function(vcf,maxSize=2000) {
  addInfoToGRanges = function(vcfinner) {
    gr = rowData(vcfinner)
    el = elementLengths(geno(vcfinner)$PR[,1])==2
    PR_NORMAL=matrix(rep(0,length(gr)*2),nc=2)
    PR_NORMAL[el] = do.call(rbind,geno(vcfinner[el])$PR[,1])
    gr$PR_NORMAL_REF=PR_NORMAL[,1]
    gr$PR_NORMAL_ALT=PR_NORMAL[,2]
    gr$PR_NORMAL_VAF = PR_NORMAL[,2]/rowSums(PR_NORMAL)
    
    el = elementLengths(geno(vcfinner)$PR[,2])==2
    PR_TUMOR=matrix(rep(0,length(gr)*2),nc=2)
    PR_TUMOR[el] = do.call(rbind,geno(vcfinner[el])$PR[,2])
    gr$PR_TUMOR_REF=PR_TUMOR[,1]
    gr$PR_TUMOR_ALT=PR_TUMOR[,2]  
    gr$PR_TUMOR_VAF = PR_TUMOR[,2]/rowSums(PR_TUMOR)
    
    el = elementLengths(geno(vcfinner)$SR[,1])==2
    SR_NORMAL=matrix(rep(0,length(gr)*2),nc=2)
    SR_NORMAL[el] = do.call(rbind,geno(vcfinner[el])$SR[,1])
    gr$SR_NORMAL_REF=SR_NORMAL[,1]
    gr$SR_NORMAL_ALT=SR_NORMAL[,2]
    gr$SR_NORMAL_VAF = SR_NORMAL[,2]/rowSums(SR_NORMAL)
    
    el = elementLengths(geno(vcfinner)$SR[,2])==2
    SR_TUMOR=matrix(rep(0,length(gr)*2),nc=2)
    SR_TUMOR[el] = do.call(rbind,geno(vcfinner[el])$SR[,2])
    gr$SR_TUMOR_REF=SR_TUMOR[,1]
    gr$SR_TUMOR_ALT=SR_TUMOR[,2]
    gr$SR_TUMOR_VAF = SR_TUMOR[,2]/rowSums(SR_TUMOR)
    
    gr$TUMOR_DP = gr$SR_TUMOR_REF + gr$SR_TUMOR_ALT + gr$PR_TUMOR_REF + gr$PR_TUMOR_ALT
    gr$NORMAL_DP = gr$SR_NORMAL_REF + gr$SR_NORMAL_ALT + gr$PR_NORMAL_REF + gr$PR_NORMAL_ALT
    
    gr$QUAL = info(vcfinner)$SOMATICSCORE
    gr$SVTYPE=info(vcfinner)$SVTYPE
    gr$SVLEN =info(vcfinner)$SVLEN
    gr$HOMLEN = info(vcfinner)$HOMLEN
    gr$HOMSEQ = info(vcfinner)$HOMSEQ
    gr$SVINSLEN = info(vcfinner)$SVINSLEN
    gr$SVINSSEQ = info(vcfinner)$SVINSSEQ
    return(gr)
  }
  end(rowData(vcf)[info(vcf)$SVTYPE=="DEL"])=info(vcf[info(vcf)$SVTYPE=="DEL"])$END
  rd = addInfoToGRanges(vcf[width(vcf)<=maxSize])
  rd2 = addInfoToGRanges(vcf[width(vcf)>maxSize])
  start(rd2)=end(rd2)
  rd3 = addInfoToGRanges(vcf[width(vcf)>maxSize])
  end(rd3)=start(rd3)
  ret = c(rd,rd2,rd3)
  names(ret)=make.unique(names(ret))
  return(ret)
}

