

### Extended Data Fgure 3a ### 

### Check links and PPIs in >7 bin
setwd('N:/Interlink_subgroup/Scripts_functions')
sapply(paste('./Functions',list.files(pattern="[.]R$", path="./Functions"),sep='/'), source);

### CW ### HEK
XLs_ori <- read.csv("N:/Interlink_subgroup/Cong_DSBSO//SS_unique_lys_crosslink.csv") ## not useful, where is the separation target-decoys?
max_score <- -log10(apply(cbind(XLs_ori$n_score_a_MS2_MS3, XLs_ori$n_score_b_MS2_MS3),1, max))
XLs_ori_CW <- cbind.data.frame(XLs_ori, max_score)
XLs_ori_CW$"ID" <- rownames(XLs_ori_CW)

### BB ### virion
XLs_ori <- read.csv("N:/Interlink_subgroup/Boris_Virion_DSSO//SS_unique_lys_crosslink.csv") ## not useful, where is the separation target-decoys?
max_score <- -log10(apply(cbind(XLs_ori$n_score_a_MS2_MS3, XLs_ori$n_score_b_MS2_MS3),1, max))
XLs_ori_BB <- cbind.data.frame(XLs_ori, max_score)
XLs_ori_BB<- GeneNameExtractXlinkX(XLs_ori_BB)
XLs_ori_BB$"ID" <- rownames(XLs_ori_BB)

### YZ ### mito 
XLs_ori <- read.csv("N:/Interlink_subgroup/Ying_DSSO//SS_unique_lys_crosslink.csv") ## not useful, where is the separation target-decoys?
max_score <- -log10(apply(cbind(XLs_ori$n_score_a_MS2_MS3, XLs_ori$n_score_b_MS2_MS3),1, max))
XLs_ori_YZ <- cbind.data.frame(XLs_ori, max_score)
XLs_ori_YZ$"ID" <- rownames(XLs_ori_YZ)



###CW###
CW_fus <- MakeInterLinkDependentGroups(XLs_ori = XLs_ori_CW, FUSED=T)
CW_sim <- FDR_simulator3(SampSize = 4860, TT = 12216, TD = 2003-69, DD=69, XLs_ori = XLs_ori_CW)
###BB###
BB_fus <- MakeInterLinkDependentGroups(XLs_ori = XLs_ori_BB, FUSED=T)
BB_sim <- FDR_simulator3(SampSize = 1318, TT = 7048-615, TD = 615, DD=0, XLs_ori = XLs_ori_BB)
###YZ###
YZ_fus <- MakeInterLinkDependentGroups(XLs_ori = XLs_ori_YZ, FUSED=T)
YZ_sim <- FDR_simulator3(SampSize = 4152, TT = 9583, TD = 681, DD=15, XLs_ori = XLs_ori_YZ)



MakeInterBarplot <- function(CW_con, titt="HEK, concatenated"){
  b <- barplot(as.matrix(CW_con),stacked=T, col=c("lightgrey", "red"), ylim=c(0,8000), 
               names.arg = c("<2 Lys per protein", "at least 2 \n Lys per protein"), ylab="inter-link count", main=titt)
  text(x = b, y=colSums(CW_con)+250,labels = c(paste("FDR = ", round(CW_con[2,1]/CW_con[1,1],2)),paste("FDR = ", round(CW_con[2,2]/CW_con[1,2],4))))
  return(as.matrix(CW_con))
}


CWse <- MakeInterBarplot(CW_fus, titt = "HEK, fused")
CWsi <- MakeInterBarplot(CW_sim, titt = "HEK, simulated")

YZse <- MakeInterBarplot(YZ_fus, titt = "mito, fused")
YZsi <- MakeInterBarplot(YZ_sim, titt = "mito, simulated")

BBse <- MakeInterBarplot(BB_fus, titt = "virion, fused")
BBsi <- MakeInterBarplot(BB_sim, titt = "virion, simulated")
