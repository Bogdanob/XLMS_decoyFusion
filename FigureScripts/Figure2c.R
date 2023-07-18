

### Figure 2



library(stringr)


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



### Cong
Links2_out_CW <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_CW, FUSED=T, Inter2=T, Grouped=T) ## context-rich
Links1_out_CW <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_CW, FUSED=T, Inter2=F, Grouped=T) ## context-poor
Old_CW <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_CW, FUSED=T, Inter2=F, Grouped=F) ## standard

FDR_LI1_CW <- MakeFDR2(FDR_cut = 0.01, sel = Links1_out_CW) ## subgroup FDR, context-rich
FDR_LI2_CW <- MakeFDR2(FDR_cut = 0.01, sel = Links2_out_CW) ## subgroup FDR, context-poor
FDR_OLD_CW <- MakeFDR2(FDR_cut = 0.01, sel = Old_CW) ## FDR, standard

Het_FDR_CW <-HarmonizeSubgroupFDR(SG1=FDR_LI1_CW, SG2 = FDR_LI2_CW, FDR_append = "inter-dependent", SG=T, XLs = XLs_ori_CW) ## combine subgroup FDRs
trad_FDR_CW <-HarmonizeSubgroupFDR(SG1=FDR_OLD_CW, FDR_append = "traditional", SG=F, XLs = XLs_ori_CW) 

ForSelsHet  <- MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_CW, removeDecoys = T) ## how many at FDR=0.01, subgrouped
ForSelsTra  <- MakeFDR3(FDR_cut = 0.01, sel = trad_FDR_CW, removeDecoys = T) ## how many at FDR=0.01, standard

### intra-dependent

CW_both <- SelectSubset(XLs_ori = XLs_ori_CW, BothIntras = T, Grouped = T) ## context-rich
CW_none <- SelectSubset(XLs_ori = XLs_ori_CW, BothIntras = F, Grouped = T) ## context-poor
FDR_CW_both <- MakeFDR2(FDR_cut = 0.01, sel=CW_both) ## subgroup FDR, context-rich 
FDR_CW_none <- MakeFDR2(FDR_cut = 0.01, sel=CW_none) ## subgroup FDR, context-poor 

Intra_FDR_CW <- HarmonizeSubgroupFDR(SG1=FDR_CW_both, SG2 = FDR_CW_none, FDR_append = "intra-dependent", SG=T, XLs = XLs_ori_CW) ## combine subgroup FDRs

FDR001_CW_inter <- length(which(MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_CW, removeDecoys = T)[[3]] < 0.01)) ## how many at FDR=0.01
FDR001_CW_trad <- length(which(MakeFDR3(FDR_cut = 0.01, sel = trad_FDR_CW, removeDecoys = T)[[3]] < 0.01)) ## how many at FDR=0.01
FDR001_CW_intra <- length(which(MakeFDR3(FDR_cut = 0.01, sel = Intra_FDR_CW, removeDecoys = T)[[3]] < 0.01)) ## how many at FDR=0.01

### Ying ### same for mito
Links2_out_YZ <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_YZ, FUSED=T, Inter2=T, Grouped=T)
Links1_out_YZ <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_YZ, FUSED=T, Inter2=F, Grouped=T)
Old_YZ <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_YZ, FUSED=T, Inter2=F, Grouped=F)

FDR_LI1_YZ <- MakeFDR2(FDR_cut = 0.01, sel = Links1_out_YZ)
FDR_LI2_YZ <- MakeFDR2(FDR_cut = 0.01, sel = Links2_out_YZ)
FDR_OLD_YZ <- MakeFDR2(FDR_cut = 0.01, sel = Old_YZ)

Het_FDR_YZ <-HarmonizeSubgroupFDR(SG1=FDR_LI1_YZ, SG2 = FDR_LI2_YZ, FDR_append = "inter-dependent", SG=T, XLs = XLs_ori_YZ)
trad_FDR_YZ <-HarmonizeSubgroupFDR(SG1=FDR_OLD_YZ, FDR_append = "traditional", SG=F, XLs = XLs_ori_YZ)

ForSelsHet  <- MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_YZ, removeDecoys = T)
ForSelsTra  <- MakeFDR3(FDR_cut = 0.01, sel = trad_FDR_YZ)



YZ_both <- SelectSubset(XLs_ori = XLs_ori_YZ, BothIntras = T, Grouped = T)
YZ_none <- SelectSubset(XLs_ori = XLs_ori_YZ, BothIntras = F, Grouped = T)
FDR_YZ_both <- MakeFDR2(FDR_cut = 0.01, sel=YZ_both)
FDR_YZ_none <- MakeFDR2(FDR_cut = 0.01, sel=YZ_none)

Intra_FDR_YZ <- HarmonizeSubgroupFDR(SG1=FDR_YZ_both, SG2 = FDR_YZ_none, FDR_append = "intra-dependent", SG=T, XLs = XLs_ori_YZ)

FDR001_YZ_inter <- length(which(MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_YZ, removeDecoys = T)[[3]] < 0.01))
FDR001_YZ_trad <- length(which(MakeFDR3(FDR_cut = 0.01, sel = trad_FDR_YZ, removeDecoys = T)[[3]] < 0.01))
FDR001_YZ_intra <- length(which(MakeFDR3(FDR_cut = 0.01, sel = Intra_FDR_YZ, removeDecoys = T)[[3]] < 0.01))

### BB ### same for virion
Links2_out_BB <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_BB, FUSED=T, Inter2=T, Grouped=T)
Links1_out_BB <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_BB, FUSED=T, Inter2=F, Grouped=T)
Old_BB <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_BB, FUSED=T, Inter2=F, Grouped=F)

FDR_LI1_BB <- MakeFDR2(FDR_cut = 0.01, sel = Links1_out_BB)
FDR_LI2_BB <- MakeFDR2(FDR_cut = 0.01, sel = Links2_out_BB)
FDR_OLD_BB <- MakeFDR2(FDR_cut = 0.01, sel = Old_BB)

Het_FDR_BB <-HarmonizeSubgroupFDR(SG1=FDR_LI1_BB, SG2 = FDR_LI2_BB, FDR_append = "inter-dependent", SG=T, XLs = XLs_ori_BB)
trad_FDR_BB <-HarmonizeSubgroupFDR(SG1=FDR_OLD_BB, FDR_append = "traditional", SG=F, XLs = XLs_ori_BB)

ForSelsHet  <- MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_BB, removeDecoys = T)
ForSelsTra  <- MakeFDR3(FDR_cut = 0.01, sel = trad_FDR_BB)



BB_both <- SelectSubset(XLs_ori = XLs_ori_BB, BothIntras = T, Grouped = T)
BB_none <- SelectSubset(XLs_ori = XLs_ori_BB, BothIntras = F, Grouped = T)
FDR_BB_both <- MakeFDR2(FDR_cut = 0.01, sel=BB_both)
FDR_BB_none <- MakeFDR2(FDR_cut = 0.01, sel=BB_none)

Intra_FDR_BB <- HarmonizeSubgroupFDR(SG1=FDR_BB_both, SG2 = FDR_BB_none, FDR_append = "intra-dependent", SG=T, XLs = XLs_ori_BB)

FDR001_BB_inter <- length(which(MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_BB, removeDecoys = T)[[3]] < 0.01))
FDR001_BB_trad <- length(which(MakeFDR3(FDR_cut = 0.01, sel = trad_FDR_BB, removeDecoys = T)[[3]] < 0.01))
FDR001_BB_intra <- length(which(MakeFDR3(FDR_cut = 0.01, sel = Intra_FDR_BB, removeDecoys = T)[[3]] < 0.01))

dissi <- c(FDR001_CW_trad, FDR001_CW_intra, FDR001_CW_inter, 
           FDR001_YZ_trad, FDR001_YZ_intra, FDR001_YZ_inter, 
           FDR001_BB_trad, FDR001_BB_intra, FDR001_BB_inter)

b <- barplot(dissi, ylim=c(0,8850),
        col=c("darkgrey", "orange", "mediumseagreen"), 
        ylab="target inter-links at FDR = 0.01", 
        names.arg = rep(c("standard", "intra-dependent", "inter-dependent"), 3),
        space=c(1,.2,.2,1,.2,.2,1,.2,.2))
grid(NA,NULL)

display <- c(FDR001_CW_trad/FDR001_CW_trad, FDR001_CW_intra/FDR001_CW_trad, FDR001_CW_inter/FDR001_CW_trad,
             FDR001_YZ_trad/FDR001_YZ_trad, FDR001_YZ_intra/FDR001_YZ_trad, FDR001_YZ_inter/FDR001_YZ_trad,
             FDR001_BB_trad/FDR001_BB_trad, FDR001_BB_intra/FDR001_BB_trad, FDR001_BB_inter/FDR001_BB_trad)
text(paste(round(display*100), "%", sep=" "), x=c(b), y=dissi+100)

