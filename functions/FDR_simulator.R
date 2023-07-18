library(viridis)
library(circlize)
library(ggplot2)
library(Rcpp)
library(ComplexHeatmap)
library(igraph)
library(stats4)


### function simulating the FPs among subgroups, dependent on several parameters ###
FDR_simulator <- function(SampSize, countIntra, countInter, probCorrectInter, AllCorrect = F, XLs_ori = XLs_ori){
  ### each Protein in the database gets a probability of intra-links
  if (AllCorrect == F){

    Aprobs <- PrepareDataFit(XLs_ori = XLs_ori, SampSize = SampSize)

    intraProbs <- sample(Aprobs[[1]], replace = F)
    wrongProbs <- sample(Aprobs[[2]], replace = F)
  }
  interProbs <- intraProbs   ### each Protein in the database gets the same probability for inter-links
  FP3 <- sample(1:SampSize,size = countIntra,prob = intraProbs, replace = T)  # assign intra-links to the proteins in the database according to the probabilities
  FP <- get_vector_with_proportion(n=countInter, p=probCorrectInter)   ### FP: make a proportional vector with 0 (incorrect) and 1 (correct) links according to p
  ## database search, correct links (1) will be placed according to the power-law derived probabilities
  ## incorrect links (0) will be placed randomly on all proteins
  FP1 <- ifelse(FP == 1, sample(SampSize, size = length(FP), replace = TRUE, prob=interProbs), 
                         sample(SampSize, size = length(FP), replace = TRUE, prob=wrongProbs))
  FP4 <- ifelse(FP == 1, sample(SampSize, size = length(FP), replace = TRUE, prob=interProbs), 
                         sample(SampSize, size = length(FP), replace = TRUE, prob=wrongProbs))
  
  Again <- (FP1 == FP4)

  ### remove spurious intra-links ###  
  while(any(Again)){
    FP1[Again] <- ifelse(FP == 1, sample(SampSize, size = length(FP), replace = TRUE, prob=interProbs), 
                  sample(SampSize, size = length(FP), replace = TRUE, prob=wrongProbs))[Again]
    FP4[Again] <- ifelse(FP == 1, sample(SampSize, size = length(FP), replace = TRUE, prob=interProbs), 
                  sample(SampSize, size = length(FP), replace = TRUE, prob=wrongProbs))[Again]
    Again <- FP1 == FP4
  }

  HaveIntra <- FP1%in%FP3 & FP4%in%FP3   # which of the inter-linked proteins have intra-links
  sorted_pasted_rows <- apply(cbind(FP1,FP4), 1, function(x) paste(sort(x), collapse=","))   # make a unique identifier for inter-links
  ## check which ones have both intra, one intra or no intra
  NoIntra <- subset(FP,(!FP1 %in% FP3 & !FP4 %in% FP3))
  BothIntra <- subset(FP,(FP1 %in% FP3 & FP4 %in% FP3))
  OneIntra <- subset(FP,(FP1 %in% FP3 | FP4 %in% FP3) &! (FP1 %in% FP3 & FP4 %in% FP3))
  OneIntraNoIntra <- c(NoIntra,OneIntra)
  # make a matrix for output
  BI_sim <- table(BothIntra==0)
  # NI_sim <- table(NoIntra==0)
  # OI_sim <- table(OneIntra==0)
  OINI_sim <- table(OneIntraNoIntra==0)
  return(matrix(ncol=2, nrow=2, c(BI_sim,OINI_sim)))
}
