

setwd('N:/Interlink_subgroup/Cong_DSBSO/')

### CW
XLs_ori <- read.csv("N:/Interlink_subgroup/Cong_DSBSO//SS_unique_lys_crosslink.csv") ## not useful, where is the separation target-decoys?
max_score <- -log10(apply(cbind(XLs_ori$n_score_a_MS2_MS3, XLs_ori$n_score_b_MS2_MS3),1, max))
XLs_ori_CW <- cbind.data.frame(XLs_ori, max_score)

### BB
XLs_ori_BB <- read.csv("N:/Interlink_subgroup/Boris_Virion_DSSO//SS_unique_lys_crosslink.csv") ## not useful, where is the separation target-decoys?
max_score <- -log10(apply(cbind(XLs_ori$n_score_a_MS2_MS3, XLs_ori$n_score_b_MS2_MS3),1, max))
XLs_ori_BB <- cbind.data.frame(XLs_ori, max_score)
XLs_ori_BB<- GeneNameExtractXlinkX(XLs_ori_BB)


## YZ
XLs_ori <- read.csv("N:/Interlink_subgroup/Ying_DSSO//SS_unique_lys_crosslink.csv") ## not useful, where is the separation target-decoys?
max_score <- -log10(apply(cbind(XLs_ori$n_score_a_MS2_MS3, XLs_ori$n_score_b_MS2_MS3),1, max))
XLs_ori_YZ <- cbind.data.frame(XLs_ori, max_score)



#### Figure 1b

pdf(pdfname <- "Figure2.pdf", width=12,height=6)

layout(matrix(ncol=6, nrow=2, c(1,8,2,8,3,4,3,5,8,6,8,7)))

### plotting
SubG_CW_sep <- MakeSubgroupFDR(XLs_ori = XLs_ori_CW, Concatenated = T)
SubG_CW_fus <- MakeSubgroupFDR(XLs_ori = XLs_ori_CW, Concatenated = F)


table(!XLs_ori_CW$gene_a == XLs_ori_CW$gene_b & XLs_ori_CW$target_decoy == "decoy")
1- 2003/14219


#SubG_CW_fus <- MakeSubgroupFDR(XLs_ori = XLs_ori_CW, Concatenated = F)
Sim_CW <- FDR_simulator(SampSize = 4860, probCorrectInter = 0.859, countIntra = 23504, countInter = 14219, XLs_ori = XLs_ori_CW) # CW

layout(matrix(ncol = 3, nrow=1, c(1:3)))
MakeSG_barplot(SubG_CW_fus = SubG_CW_sep, titt="concatenated search", ylim=c(0,16000), Simulated = F)
MakeSG_barplot(SubG_CW_fus = Sim_CW, titt="simulated search", ylim=c(0,16000), Simulated = T)
MakeSG_barplot(SubG_CW_fus = SubG_CW_fus,titt = "fused search", Simulated = F, ylim=c(0,16000))



### Extendend Data Figure 2b
#### make a heatmap comparing intra-link count and inter-link prior probability
Returnis <- c(1:30)*1000
Returnis2 <- c(1:10)*1000
result_matrix <- matrix(NA, nrow=length(Returnis), ncol=length(Returnis2))
for (i in 1:length(Returnis)) {
  for (j in 1:length(Returnis2)) {
    temp<- FDR_simulator(SampSize = Returnis2[j], countInter = 12219, probCorrectInter = 0.804, countIntra = Returnis[i], XLs_ori = XLs_ori_CW)
    result_matrix[i, j] <- temp[2,1]/(temp[1,1]+temp[2,1])
  }
  cat(Returnis[i], "\r")
}

rownames(result_matrix) <- Returnis/100
colnames(result_matrix) <- Returnis2/100

image(x=Returnis/100, y=Returnis2/100,z=result_matrix, col = viridis(20), 
      xlab="intra-link count x 100", ylab="database entries x 100", 
      main="#inter-links: 12,219, prob. for correct inter-link = 0.8")
### for Legend
col_fun = colorRamp2(seq(from=min(result_matrix), to = max(result_matrix), b=0.001)
                     , viridis(length(seq(from=min(result_matrix), to = max(result_matrix), b=0.001))))

#seq(from=min(result_matrix), to = max(result_matrix), b=0.001)
lgd = Legend(col_fun = col_fun, title = "FP [%]")
pushViewport(viewport(width = 1, height = 1))
draw(lgd, x = unit(0.6, "cm"), y = unit(0.6, "cm"), just = c("left", "bottom"))

### Extendend Data Figure 2c

SubG_YZ_sep <- MakeSubgroupFDR(XLs_ori = XLs_ori_YZ, Concatenated = T)
SubG_YZ_fus <- MakeSubgroupFDR(XLs_ori = XLs_ori_YZ, Concatenated = F)
Sim_YZ <- FDR_simulator(SampSize = 4152, probCorrectInter =1-(696/(9583-696)), countIntra = 6250, countInter = 9583-696, XLs_ori = XLs_ori_YZ) # YZ

SubG_BB_sep <- MakeSubgroupFDR(XLs_ori = XLs_ori_BB, Concatenated = T)
SubG_BB_fus <- MakeSubgroupFDR(XLs_ori = XLs_ori_BB, Concatenated = F)
Sim_BB <- FDR_simulator(SampSize = 2500, probCorrectInter =0.859, countIntra = 9313, countInter = 3633, XLs_ori = XLs_ori_BB) # BB

FDRs_sep <- c(SubG_CW_sep[2,1]/SubG_CW_sep[1,1], SubG_BB_sep[2,1]/SubG_BB_sep[1,1], SubG_YZ_sep[2,1]/SubG_YZ_sep[1,1])
FDRs_fus <- c(SubG_CW_fus[2,1]/SubG_CW_fus[1,1], SubG_BB_fus[2,1]/SubG_BB_fus[1,1], SubG_YZ_fus[2,1]/SubG_YZ_fus[1,1])
FPRs <- c(Sim_CW[2,1]/Sim_CW[1,1], Sim_BB[2,1]/Sim_BB[1,1], Sim_YZ[2,1]/Sim_YZ[1,1])

FDRs_sep_2 <- c(SubG_CW_sep[2,2]/SubG_CW_sep[1,2], SubG_BB_sep[2,2]/SubG_BB_sep[1,2], SubG_YZ_sep[2,2]/SubG_YZ_sep[1,2])
FDRs_fus_2 <- c(SubG_CW_fus[2,2]/SubG_CW_fus[1,2], SubG_BB_fus[2,2]/SubG_BB_fus[1,2], SubG_YZ_fus[2,2]/SubG_YZ_fus[1,2])
FPRs_2 <- c(Sim_CW[2,2]/Sim_CW[1,2], Sim_BB[2,2]/Sim_BB[1,2], Sim_YZ[2,2]/Sim_YZ[1,2])


plot(cbind(FDRs_fus,FPRs), xlim=c(0,0.2), ylim=c(0,0.2),
     xlab="FDR in subgroup", ylab="FPR, simulated", frame=F, col="orange", pch=20, cex=2)
points(cbind(FDRs_sep,FPRs), xlim=c(0,0.2), ylim=c(0,0.2),
     xlab="FDR in subgroup", ylab="FPR, simulated", frame=F, col="mediumseagreen", pch=20, cex=2)
abline(a=0, b=1)




graphics.off()
system(paste("open", pdfname))




