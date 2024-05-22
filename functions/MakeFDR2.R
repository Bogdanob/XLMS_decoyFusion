MakeFDR2 <-function(FDR_cut, sel){
#    FDR_cut <- 0.02
 #   sel=MaxBI
 # FDR_cut = 0.01; sel = Links1_out_CW
  
  #  FDR_cut <- 0.05
  #  sel <-Het_FDR_MXR
  
 # sel=MaxBIc; FDR_cut = 0.01
  
  Revs <- sel$target_decoy == "decoy" ## target-decoys + decoy-decoy
  DRevs <- grepl(x = sel$protein_a_cln, pattern="###RND###") & grepl(x = sel$protein_b_cln, pattern="###RND###")  ## decoy-decoy
  Revs <- Revs &! DRevs ## target-decoys
  
  
  
  All_scores <- sel$max_score
  FW_scores <- sel$max_score[!Revs]
  Rev_scores <- sel$max_score[Revs] ## &!DRevs
  Rev2_scores<- sel$max_score[DRevs]
  
  names(All_scores) <- sel$ID
  
  
  SAll <- sort(All_scores, decreasing = T)
  SR <- sort(Rev_scores, decreasing = T)
  SF <- sort(FW_scores, decreasing = T)
  SR2 <- sort(Rev2_scores, decreasing = T)
  
  FDR<-NULL
  for (i in 1:length(SAll)){ ## run through all IDs
    if(any(SR>=SAll[i])){
      #      i=3000
      dec <- max(which(SR>SAll[i])) ## target-decoys
      
      dec2 <-max(which(SR2>SAll[i])) ## decoy-decoy
      
      hit <- max(which(SF>=SAll[i])) ## target-target
      
      if (dec2 != "-Inf"){
        decOverhit <- (dec-dec2)/hit  
      }
      if (dec2 == "-Inf"){
        decOverhit <- (dec)/hit  
      }
      
      FDR[i]<-decOverhit
      
    }
    cat(i/length(SAll)*100,"\r")
    
  }
  #FDR_cut <- 0.01
  Recalled <- c(1:length(FDR))
  Cut<- SAll[max(which(FDR<=FDR_cut))]
  if(is.na(Cut)){
    Cut <- 0
  }
  FDR[FDR=="-Inf"] <- NA
  
  if (min(FDR, na.rm=T) == "Inf"){
    FDR <- rep(0,length(SAll))
  }
  FDR[is.na(FDR)] <- min(FDR, na.rm=T)

  # abline(v = Cut)
  return(list(Cut,Recalled,FDR, names(SAll)))
  # plot(y=Recalled, x=FDR, type="l", main=tile)
}
