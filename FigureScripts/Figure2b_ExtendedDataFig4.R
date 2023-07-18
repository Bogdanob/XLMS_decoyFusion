

library(stringr)

### inter-link dependent grouping ###


### Check links and PPIs in >7 bin
setwd('N:/Interlink_subgroup/Scripts_functions')
sapply(paste('./Functions',list.files(pattern="[.]R$", path="./Functions"),sep='/'), source);

### CW
XLs_ori <- read.csv("N:/Interlink_subgroup/Cong_DSBSO//SS_unique_lys_crosslink.csv") ## not useful, where is the separation target-decoys?
max_score <- -log10(apply(cbind(XLs_ori$n_score_a_MS2_MS3, XLs_ori$n_score_b_MS2_MS3),1, max))
XLs_ori_CW <- cbind.data.frame(XLs_ori, max_score)
XLs_ori_CW$"ID" <- rownames(XLs_ori_CW)

### BB
XLs_ori <- read.csv("N:/Interlink_subgroup/Boris_Virion_DSSO//SS_unique_lys_crosslink.csv") ## not useful, where is the separation target-decoys?
max_score <- -log10(apply(cbind(XLs_ori$n_score_a_MS2_MS3, XLs_ori$n_score_b_MS2_MS3),1, max))
XLs_ori_BB <- cbind.data.frame(XLs_ori, max_score)
XLs_ori_BB<- GeneNameExtractXlinkX(XLs_ori_BB)
XLs_ori_BB$"ID" <- rownames(XLs_ori_BB)

## YZ
XLs_ori <- read.csv("N:/Interlink_subgroup/Ying_DSSO//SS_unique_lys_crosslink.csv") ## not useful, where is the separation target-decoys?
max_score <- -log10(apply(cbind(XLs_ori$n_score_a_MS2_MS3, XLs_ori$n_score_b_MS2_MS3),1, max))
XLs_ori_YZ <- cbind.data.frame(XLs_ori, max_score)
XLs_ori_YZ$"ID" <- rownames(XLs_ori_YZ)


###CW###
CW_con <- MakeInterLinkDependentGroups(XLs_ori = XLs_ori_CW, FUSED=F)
CW_fus <- MakeInterLinkDependentGroups(XLs_ori = XLs_ori_CW, FUSED=T)
###BB###
BB_con <- MakeInterLinkDependentGroups(XLs_ori = XLs_ori_BB, FUSED=F)
BB_fus <- MakeInterLinkDependentGroups(XLs_ori = XLs_ori_BB, FUSED=T)
###YZ###
YZ_con <- MakeInterLinkDependentGroups(XLs_ori = XLs_ori_YZ, FUSED=F)
YZ_fus <- MakeInterLinkDependentGroups(XLs_ori = XLs_ori_YZ, FUSED=T)






pdf(pdfname<- "Figure2b_ExtDataFig4.pdf", width=12, height=4)
layout(matrix(ncol=3, nrow=1, c(1:3)))


### cw 
Links2_out_CW <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_CW, FUSED=T, Inter2=T, Grouped=T)
Links1_out_CW <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_CW, FUSED=T, Inter2=F, Grouped=T)
Old_CW <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_CW, FUSED=T, Inter2=F, Grouped=F)

FDR_LI1_CW <- MakeFDR2(FDR_cut = 0.01, sel = Links1_out_CW)
FDR_LI2_CW <- MakeFDR2(FDR_cut = 0.01, sel = Links2_out_CW)
FDR_OLD_CW <- MakeFDR2(FDR_cut = 0.01, sel = Old_CW)

Het_FDR_CW <-HarmonizeSubgroupFDR(SG1=FDR_LI1_CW, SG2 = FDR_LI2_CW, FDR_append = "inter-dependent", SG=T, XLs = XLs_ori_CW)
trad_FDR_CW <-HarmonizeSubgroupFDR(SG1=FDR_OLD_CW, FDR_append = "traditional", SG=F, XLs = XLs_ori_CW)

### intra-dependent 
CW_both <- SelectSubset(XLs_ori = XLs_ori_CW, BothIntras = T, Grouped = T)
CW_none <- SelectSubset(XLs_ori = XLs_ori_CW, BothIntras = F, Grouped = T)

FDR_CW_both <- MakeFDR2(FDR_cut = 0.01, sel=CW_both)
FDR_CW_none <- MakeFDR2(FDR_cut = 0.01, sel=CW_none)

Intra_FDR_CW <- HarmonizeSubgroupFDR(SG1=FDR_CW_both, SG2 = FDR_CW_none, FDR_append = "intra-dependent", SG=T, XLs = XLs_ori_CW)
Intra_FDR_CW_filtered  <- HarmonizeSubgroupFDR(SG1=FDR_CW_both, SG2 = NULL, FDR_append = "intra-dependent", SG=T, XLs = XLs_ori_CW)



### CW
plot(y=MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_CW, removeDecoys = T)[[2]], x=MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_CW)[[3]], type="l", 
     xlim=c(0,0.1),xlab="FDR, \n (TD-DD)/TT", ylab="Number of recalled target inter-links", col="mediumseagreen", lwd=2, main="HEK")

lines(y=MakeFDR3(FDR_cut = 0.01, sel=Intra_FDR_CW, removeDecoys = T)[[2]], x= MakeFDR3(FDR_cut=0.01, sel=Intra_FDR_CW)[[3]], col="orange", lwd=2)
lines(y=MakeFDR3(FDR_cut = 0.01, sel = trad_FDR_CW)[[2]], x=MakeFDR3(FDR_cut = 0.01, sel = trad_FDR_CW)[[3]], type="l", col="darkgrey", lwd=2)
grid(NULL,NULL)
abline(v=0.01, lty=3, col="darkgrey")
legend("bottomright", legend = c("standard", "inter-dependent groups", "intra-dependent groups"), col=c("darkgrey", "mediumseagreen", "orange"), lwd=3, cex=1)


### yz ### mito
Links2_out_YZ <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_YZ, FUSED=T, Inter2=T, Grouped=T)
Links1_out_YZ <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_YZ, FUSED=T, Inter2=F, Grouped=T)
Old_YZ <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_YZ, FUSED=T, Inter2=F, Grouped=F)

FDR_LI1_YZ <- MakeFDR2(FDR_cut = 0.01, sel = Links1_out_YZ)
FDR_LI2_YZ <- MakeFDR2(FDR_cut = 0.01, sel = Links2_out_YZ)
FDR_OLD_YZ <- MakeFDR2(FDR_cut = 0.01, sel = Old_YZ)

Het_FDR_YZ <-HarmonizeSubgroupFDR(SG1=FDR_LI1_YZ, SG2 = FDR_LI2_YZ, FDR_append = "inter-dependent", SG=T, XLs = XLs_ori_YZ)
trad_FDR_YZ <-HarmonizeSubgroupFDR(SG1=FDR_OLD_YZ, FDR_append = "traditional", SG=F, XLs = XLs_ori_YZ)

### intra-dependent 
YZ_both <- SelectSubset(XLs_ori = XLs_ori_YZ, BothIntras = T, Grouped = T)
YZ_none <- SelectSubset(XLs_ori = XLs_ori_YZ, BothIntras = F, Grouped = T)
FDR_YZ_both <- MakeFDR2(FDR_cut = 0.01, sel=YZ_both)
FDR_YZ_none <- MakeFDR2(FDR_cut = 0.01, sel=YZ_none)
Intra_FDR_YZ <- HarmonizeSubgroupFDR(SG1=FDR_YZ_both, SG2 = FDR_YZ_none, FDR_append = "intra-dependent", SG=T, XLs = XLs_ori_YZ)

### YZ
plot(y=MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_YZ, removeDecoys = T)[[2]], x=MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_YZ)[[3]], type="l", 
     xlim=c(0,0.1),xlab="FDR, \n (TD-DD)/TT", ylab="Number of recalled target inter-links", col="mediumseagreen", lwd=2, main="mito")
lines(y=MakeFDR3(FDR_cut = 0.01, sel=Intra_FDR_YZ, removeDecoys = T)[[2]], x= MakeFDR3(FDR_cut=0.01, sel=Intra_FDR_YZ)[[3]], col="orange", lwd=2)
lines(y=MakeFDR3(FDR_cut = 0.01, sel = trad_FDR_YZ)[[2]], x=MakeFDR3(FDR_cut = 0.01, sel = trad_FDR_YZ)[[3]], type="l", col="darkgrey", lwd=2)
grid(NULL,NULL)
abline(v=0.01, lty=3, col="darkgrey")
legend("bottomright", legend = c("standard", "inter-dependent groups", "intra-dependent groups"), col=c("darkgrey", "mediumseagreen", "orange"), lwd=3, cex=1)

### BB ### virion
Links2_out_BB <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_BB, FUSED=T, Inter2=T, Grouped=T)
Links1_out_BB <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_BB, FUSED=T, Inter2=F, Grouped=T)
Old_BB <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_BB, FUSED=T, Inter2=F, Grouped=F)

FDR_LI1_BB <- MakeFDR2(FDR_cut = 0.01, sel = Links1_out_BB)
FDR_LI2_BB <- MakeFDR2(FDR_cut = 0.01, sel = Links2_out_BB)
FDR_OLD_BB <- MakeFDR2(FDR_cut = 0.01, sel = Old_BB)

Het_FDR_BB <-HarmonizeSubgroupFDR(SG1=FDR_LI1_BB, SG2 = FDR_LI2_BB, FDR_append = "inter-dependent", SG=T, XLs = XLs_ori_BB)
trad_FDR_BB <-HarmonizeSubgroupFDR(SG1=FDR_OLD_BB, FDR_append = "traditional", SG=F, XLs = XLs_ori_BB)

### intra-dependent 
BB_both <- SelectSubset(XLs_ori = XLs_ori_BB, BothIntras = T, Grouped = T)
BB_none <- SelectSubset(XLs_ori = XLs_ori_BB, BothIntras = F, Grouped = T)
FDR_BB_both <- MakeFDR2(FDR_cut = 0.01, sel=BB_both)
FDR_BB_none <- MakeFDR2(FDR_cut = 0.01, sel=BB_none)
Intra_FDR_BB <- HarmonizeSubgroupFDR(SG1=FDR_BB_both, SG2 = FDR_BB_none, FDR_append = "intra-dependent", SG=T, XLs = XLs_ori_BB)



### BB
plot(y=MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_BB, removeDecoys = T)[[2]], x=MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_BB)[[3]], type="l", 
     xlim=c(0,0.1),xlab="FDR, \n (TD-DD)/TT", ylab="Number of recalled target inter-links", col="mediumseagreen", lwd=2, main="virion")
lines(y=MakeFDR3(FDR_cut = 0.01, sel=Intra_FDR_BB, removeDecoys = T)[[2]], x= MakeFDR3(FDR_cut=0.01, sel=Intra_FDR_BB)[[3]], col="orange", lwd=2)
lines(y=MakeFDR3(FDR_cut = 0.01, sel = trad_FDR_BB)[[2]], x=MakeFDR3(FDR_cut = 0.01, sel = trad_FDR_BB)[[3]], type="l", col="darkgrey", lwd=2)
grid(NULL,NULL)
abline(v=0.01, lty=3, col="darkgrey")
legend("bottomright", legend = c("standard", "inter-dependent groups", "intra-dependent groups"), col=c("darkgrey", "mediumseagreen", "orange"), lwd=3, cex=1)

graphics.off()
system(paste("open", pdfname))

