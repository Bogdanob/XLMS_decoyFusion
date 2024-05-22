

FDR_simulator3 <- function(SampSize, TT, TD, DD, XLs_ori=XLs_ori){
  ### each Protein in the database gets a probability of intra-links
#  XLs_ori =XLs_ori_CW
#    SampSize = 4860; TT = 12218; TD = 2003-69; DD=69
  Aprobs <- PrepareDataFit(XLs_ori = XLs_ori, SampSize = SampSize)
  intraProbs <- sample(Aprobs[[1]], replace = F)
  wrongProbs <- sample(Aprobs[[2]], replace = F)
  
  interProbs <- sample(intraProbs, replace = F) ## randomize probabilities

  vec <- c(rep("0_0",DD),rep("1_0", TD), rep("1_1",TT))
  
  vec <- sample(vec, replace=F)
  
  
  FP_1 <- unlist(lapply(str_split(vec, "_"), "[[",1))   
  FP_2 <- unlist(lapply(str_split(vec, "_"), "[[",2))
  
  ### FP: make a proportional vector with 0 (incorrect) and 1 (correct) links according to p
  ## database search, correct links (1) will be placed according to the power-law derived probabilities
  ## incorrect links (0) will be placed randomly on all proteins

  FP_firstPos <- ifelse(FP_1 == 1, sample(SampSize, size = length(FP_1), replace = TRUE, prob=interProbs), 
                                   sample(SampSize, size = length(FP_1), replace = TRUE, prob=wrongProbs))
  tempMod   <- PrepareDataFit2(XLs_ori=XLs_ori, SampSize=SampSize)[[3]]
  coeffi_zipf <- PrepareDataFit2(XLs_ori=XLs_ori, SampSize=SampSize)[[4]]
  FP_firstPos_table <- sort(table(FP_firstPos), decreasing=T)
  
  
  FP_secondPos <- NULL
  for (i in seq_along(FP_firstPos_table)){
    set.seed(i)
    nInts <- ceiling(2^(tempMod$coefficients[2]*log2(FP_firstPos_table[i])+tempMod$coefficients[1]))+5
    corrProb <- exp(lzipf(coeffi_zipf,  nInts))
    corrProb <- c(corrProb, rep(0, SampSize - length(corrProb)))
    corrProb <- sample(corrProb, replace = F)
    
    sel <- FP_firstPos %in% names(FP_firstPos_table[i])
    FP_secondPos[sel] <- ifelse(FP_2[sel] == 1, sample(SampSize, size = length(FP_2[sel]), replace = TRUE, prob=corrProb), 
                                                sample(SampSize, size = length(FP_2[sel]), replace = TRUE, prob=wrongProbs))
    cat(i, "\r")
  }
  
  sorted_pasted_rows <- apply(cbind(FP_firstPos,FP_secondPos), 1, function(x) paste(sort(x), collapse=","))   # make a unique identifier for inter-links
  FP_real <- sorted_pasted_rows
  
  ## check how many PPIs with FPs
  PPI1 <- names(table(FP_real)[table(FP_real)== 1])
  PPI2 <- names(table(FP_real)[table(FP_real)>1])
  PPI3 <- names(table(FP_real)[table(FP_real)>2])
  PPI4 <- names(table(FP_real)[table(FP_real)>3])
  
  TPs1 <- length(FP_real[FP_1==1 & FP_2==1][FP_real[FP_1==1 & FP_2==1] %in% PPI1])
  FPs1 <- length(FP_real[!(FP_1==1 & FP_2==1)][FP_real[!(FP_1==1 & FP_2==1)] %in% PPI1])
  TPs2 <- length(FP_real[FP_1==1 & FP_2==1][FP_real[FP_1==1 & FP_2==1] %in% PPI2])
  FPs2 <- length(FP_real[!(FP_1==1 & FP_2==1)][FP_real[!(FP_1==1 & FP_2==1)] %in% PPI2])
  return(matrix(ncol=2, nrow=2, c(TPs1,FPs1, TPs2,FPs2)))
}
