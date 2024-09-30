### Makes FDR based on PEPs, removeDecoys: remove decoys from the output?, NO_Subgroups, subgrouped PEP data or non-subgrouped data?

MakePEP_FDR3 <- function(FDR_cut, sel, removeDecoys = T, NO_SUBGROUPS = F, MSAnnika=F){
  
  if(!MSAnnika){
    Revs <- sel$target_decoy == "decoy" ## target-decoys + decoy-decoy
    DRevs <- grepl(x = sel$protein_a_cln, pattern="###RND###") & grepl(x = sel$protein_b_cln, pattern="###RND###")  ## decoy-decoy
    Revs <- Revs &! DRevs ## target-decoys  
  }
  
  if(MSAnnika){
    Revs <- sel$targetdecoy == T ## target-decoys
    DRevs <- sel$decoydecoy == T ## decoy-decoy
  }
  
  
  if (NO_SUBGROUPS){
      ColSel <- which(grepl("Combined.Score", names(sel)) | grepl("max_score", names(sel)))
      All_scores <- sel[,ColSel]
      FW_scores <- sel[!Revs&!DRevs,ColSel]
      Rev_scores <- sel[Revs,ColSel]
      Rev2_scores<- sel[DRevs,ColSel]
      
      names(All_scores) <- sel$ID
  }
  
  if (NO_SUBGROUPS == F){
    ColSel <- which(grepl("PEP",names(sel)))
    All_scores <- -log10(sel[,ColSel])
    FW_scores <- -log10(sel[!Revs&!DRevs,ColSel])
    Rev_scores <- -log10(sel[Revs,ColSel])
    Rev2_scores<- -log10(sel[DRevs,ColSel])
    
    names(All_scores) <- sel$ID
    
  }
  
  SAll <- sort(All_scores, decreasing = T)
  SR <- sort(Rev_scores, decreasing = T)
  SF <- sort(FW_scores, decreasing = T)
  SR2 <- sort(Rev2_scores, decreasing = T)
  
  FDR<-NULL
  for (i in 1:length(SAll)){ ## run through all IDs
    if(any(SR>=SAll[i])){
   #  i=3000
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
  FDR[is.na(FDR)] <- min(FDR, na.rm=T)
  
  if (removeDecoys == T){
    FDR <- FDR[!(names(SAll) %in% sel$ID[sel$target_decoy =="decoy"])]
    Recalled <- c(1:length(FDR))
    SAll <- SAll[!(names(SAll) %in% sel$ID[sel$target_decoy =="decoy"])]
 # cat("here you go")
  }
  
  
  Out1 <- cbind.data.frame(names(SAll), FDR)
  names(Out1) <- c("ID", "FDR")
  Out2 <- merge(Out1, sel, by="ID", all=F)
  
  return(list(Cut,Recalled,FDR, names(SAll), Out2))
}


