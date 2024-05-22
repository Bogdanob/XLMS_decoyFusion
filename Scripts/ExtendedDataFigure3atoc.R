
### Extended Data Figure3a-c

setwd('N:/Interlink_subgroup/Scripts_functions')
sapply(paste('./Functions',list.files(pattern="[.]R$", path="./Functions"),sep='/'), source);

setwd('N:/Interlink_subgroup/Cong_DSBSO/')


XLs_ori <- read.csv("N:/Interlink_subgroup/Boris_Virion_DSSO//SS_unique_lys_crosslink.csv") ## not useful, where is the separation target-decoys?
max_score <- -log10(apply(cbind(XLs_ori$n_score_a_MS2_MS3, XLs_ori$n_score_b_MS2_MS3),1, max))
XLs_ori_CW <- cbind.data.frame(XLs_ori, max_score)
XLs_ori_CW<- GeneNameExtractXlinkX(XLs_ori_CW)



SubG_CW_sep <- MakeSubgroupFDR(XLs_ori = XLs_ori_CW, Concatenated = T)

Sim_CW <- FDR_simulator(SampSize = 1318, probCorrectInter =1-(615/(7043-615)), countIntra = 6513, countInter = 7043, XLs_ori = XLs_ori_CW) # CW
SubG_CW_fus <- MakeSubgroupFDR(XLs_ori = XLs_ori_CW, Concatenated = F)


FDR_simulator3(SampSize = 1318, TT = 7048-615, TD = 615, DD=0, XLs_ori = XLs_ori_BB)


CW_con <- MakeInterLinkDependentGroups(XLs_ori = XLs_ori_CW, FUSED=F)
CW_fus <- MakeInterLinkDependentGroups(XLs_ori = XLs_ori_CW, FUSED=T)
CW_sim <- FDR_simulator3(SampSize = 1318, TT = 7048-615, TD = 615, DD=0, XLs_ori = XLs_ori_CW)



layout(matrix(ncol=3, nrow=1, c(1:3)))

### concatenated

toPlot_con <- cbind(as.matrix(rowSums(SubG_CW_sep)), SubG_CW_sep[,c(2,1)], CW_con)

barplot(toPlot_con, main="concatenated database", ylim=c(0,10000),beside=F,
        col=c("lightgrey", "red"), 
        ylab="inter-link count", space =c(0.2,1,0.2,1,0.2,0.5))

legend(legend = c("target", "decoy"), "topright", col = c("lightgrey", "red"), pch = 15, cex=1.5)


### simulated
toPlot_sim <- cbind(Sim_CW[,c(2,1)], CW_sim)

barplot(toPlot_sim, main="simulation", ylim=c(0,10000),beside=F,
        col=c("lightgrey", "red"), 
        ylab="inter-link count", space =c(0.2,0.2,1,0.2,0.2))

legend(legend = c("true positives", "false positives"), "topright", col = c("lightgrey", "red"), pch = 15, cex=1.5)

### fused
toPlot_fus <- cbind(as.matrix(rowSums(SubG_CW_fus)), SubG_CW_fus[,c(2,1)], CW_fus)

barplot(toPlot_fus, main="fused database", ylim=c(0,10000),beside=F,
        col=c("lightgrey", "red"), 
        ylab="inter-link count", space =c(0.2,1,0.2,1,0.2,0.5))


legend(legend = c("target", "decoy"), "topright", col = c("lightgrey", "red"), pch = 15, cex=1.5)
