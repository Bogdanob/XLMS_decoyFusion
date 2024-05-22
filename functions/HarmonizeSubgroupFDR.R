

### 
HarmonizeSubgroupFDR <- function(SG1= FDR_LI1_MXR,SG2= FDR_LI2_MXR, FDR_append = "inter-dependent", SG=T, XLs=XLs_ori_P, Combined=F){
  
  
#  SG1=FDR_LI1_CW; SG2 = FDR_LI2_CW; FDR_append = "inter-dependent"; SG=T; XLs = XLs_ori_P
#  SG1=FDR_OLD_MXR 
#  SG2 = FDR_LI2_MXR 
#  FDR_append = "trad" 
#  XLs = XLs_ori_P 
#  SG = F
#  SG1=Cbd_MXR; FDR_append = "traditional_conc"; SG=F; XLs = XLs_ori_P; Combined = T
  
  if(SG==T){
    names(SG1[[3]]) <- SG1[[4]]
    names(SG2[[3]]) <- SG2[[4]]
    All <- sort(c(unlist(SG1[[3]]), unlist(SG2[[3]])), decreasing=F)    
  }  
  
  if(SG==F){
    names(SG1[[3]]) <- SG1[[4]]
    All <- sort(unlist(SG1[[3]]), decreasing=F)    
  }
  
  if (Combined == F){
    ForOut <- XLs[(XLs$protein_a_cln != XLs$protein_b_cln),]
    All <- cbind.data.frame(All, names(All))
    names(All) <- c(paste("FDR",FDR_append, sep="_"), "ID")
    ForOut$ID <- rownames(ForOut)
    Output <- merge(ForOut,All, by="ID")
    }
  
  if (Combined == T){
    ForOut <- XLs
    All <- cbind.data.frame(All, names(All))
    names(All) <- c(paste("FDR",FDR_append, sep="_"), "ID")
    ForOut$ID <- rownames(ForOut)
    Output <- merge(ForOut,All, by="ID")
  }
  
  return(Output)
}

