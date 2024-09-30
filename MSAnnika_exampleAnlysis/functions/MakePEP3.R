### claculate PEPs for MSAnnika, sel = input data, MSAnnika input = T/F

MakePEP3 <- function(sel, MSAnnika=F){
  
  if(!MSAnnika){
  Revs <- sel$target_decoy == "decoy" # Identify decoys
  }
  
  if(MSAnnika){
    Revs <- sel$targetdecoy == T
  }
  
  FW_scores <- sel$max_score[!Revs] # extract targets
  Rev_scores <- sel$max_score[Revs] # extract decoys
  
  est <- approxfun(density(FW_scores, bw="nrd0",adjust=2, from=min(sel$max_score), to=max(sel$max_score))) # gaussian kernel densities, targets
  
  if (length(Rev_scores) != 0){
    rest <- approxfun(density(Rev_scores, bw="nrd0",adjust=2, from=min(sel$max_score), to=max(sel$max_score))) # gaussian kernel densities, decoys  
  }
  
  if (length(Rev_scores) == 0){
    rest <- function(x) {return(0)} 
  }
  
  
  targets <- length(FW_scores) # number of targets
  decoys <- length(Rev_scores) # number of decoys
  total <- targets + decoys # all
  
  Pdecoy <- decoys / total # prior P decoy
  Ptarget <- targets / total # prior P target
  
  PEPs <- function(x) {   
    P_score_decoy = rest(x)
    P_score_target = est(x)
    numerator = P_score_decoy * Pdecoy
    denominator = P_score_decoy * Pdecoy + P_score_target * Ptarget
    return(numerator / denominator)
  }
  sel$PEP <- PEPs(sel$max_score)

  return(sel)
}


#plot(sel$max_score, sel$PEP)
