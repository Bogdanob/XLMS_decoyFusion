

### new Figure 1 ###



#### Figure 1 plots and analyses ####

### Max ### ground-truth comparison
setwd('N:/Interlink_subgroup/Scripts_functions')
sapply(paste('./Functions',list.files(pattern="[.]R$", path="./Functions"),sep='/'), source);

setwd("N:/Interlink_subgroup/Max_plates")

Truth <- read.csv(file = "all_plates_rerun_selectedfiles_Claasen_etAl_XlinkX.csv")
dim(Truth)
### how to deal with contaminants # remove
Truth <- Truth[!(grepl(Truth$protein_a, pattern="contam") | grepl(Truth$protein_b, pattern="contam")),] #& !(Truth$group_a != "" | Truth$group_b != "")),]

Truth <- Truth[!Truth$Truth == "CONTAMINANT",]

XLs_ori_1 <- Truth 
max_score <- -log10(apply(cbind(XLs_ori_1$n_score_a_MS2_MS3, XLs_ori_1$n_score_b_MS2_MS3),1, max))
XLs_ori_1 <- cbind.data.frame(XLs_ori_1, max_score)

XLs_ori_P <- XLs_ori_1
XLs_ori_P$ID <- rownames(XLs_ori_P)


Old_MXR <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_P, FUSED=F, Inter2=F, Grouped=F)


### make subgroup FDR ### inter-link dependent ### fused
Old_FDR_MXR <- MakeFDR2(sel=Old_MXR, FDR_cut = 0.01)
Cbd_MXR <- MakeFDR2(sel = XLs_ori_P, FDR_cut = 0.01)

### intra-link dependent Groups ### fused search
MaxOld <- SelectSubset(XLs_ori = XLs_ori_P, BothIntras = F, Grouped = F, Concatenated = F)

### intra-link dependent groups ### concatenated search 
MaxOldc <- SelectSubset(XLs_ori = XLs_ori_P, BothIntras = F, Grouped = F, Concatenated = T)




### make subgroup FDR ### intra-dependent ### fused
FDR_MxR_old      <- MakeFDR2(sel=MaxOld, FDR_cut = 0.01)

### make subgroup FDR ### intra-dependent ### concatenated
FDR_MxR_oldc      <- MakeFDR2(sel=MaxOldc, FDR_cut = 0.01)

### inter-link dependent groups ### fused search 
trad_FDR_MXR <-HarmonizeSubgroupFDR(SG1=Old_FDR_MXR, FDR_append = "traditional", SG=F, XLs = XLs_ori_P)
#Cbd_FDR_MXR <-HarmonizeSubgroupFDR(SG1=Cbd_MXR, FDR_append = "combined", SG=F, XLs = XLs_ori_P)



### inter-link dependent groups ### concatenated search 
#Het_FDR_MXRc  <-HarmonizeSubgroupFDR(SG1=FDR_LI1_MXRc, SG2 = FDR_LI2_MXRc, FDR_append = "inter-dependent_conc", SG=T, XLs = XLs_ori_P)
trad_FDR_MXRc <-HarmonizeSubgroupFDR(SG1=Old_FDR_MXR, FDR_append = "traditional_conc", SG=F, XLs = XLs_ori_P)
trad_FDR_MXRc <- trad_FDR_MXRc[!(trad_FDR_MXRc$Truth =="DECOY" | trad_FDR_MXRc$Truth =="INTRALINK"),]

Cbd_FDR_MXRc <-HarmonizeSubgroupFDR(SG1=Cbd_MXR, FDR_append = "traditional_combined", SG=F, XLs = XLs_ori_P, Combined = T)
Cbd_FDR_MXRc <- Cbd_FDR_MXRc[!(Cbd_FDR_MXRc$Truth =="DECOY" | Cbd_FDR_MXRc$Truth =="INTRALINK"),]


spec <- NULL
for (i in 1:5){
  i <- i/100
  cat(i)
  a <- sum(trad_FDR_MXRc$FDR_traditional_conc < i)
  b <- sum(trad_FDR_MXRc$FDR_traditional_conc < i & trad_FDR_MXRc$Truth == F)

  spec[i*100] <- 1-(b/a)
  names(spec)[i*100] <- i*100
  }

sens <- NULL
for (i in 1:5){
  i <- i/100
  cat(i)
  a <- sum(trad_FDR_MXRc$FDR_traditional_conc < i & trad_FDR_MXRc$Truth == T)
  b <- sum(trad_FDR_MXRc$Truth == T)
  
  sens[i*100] <- (b-a)/b
  names(sens)[i*100] <- i*100
}


spec_Cbd <- NULL
for (i in 1:5){
  i <- i/100
  cat(i)
  a <- sum(Cbd_FDR_MXRc$FDR_traditional_combined < i)
  b <- sum(Cbd_FDR_MXRc$FDR_traditional_combined < i & Cbd_FDR_MXRc$Truth == F)
  
  spec_Cbd[i*100] <- 1-(b/a)
  names(spec_Cbd)[i*100] <- i*100
}

sens_Cbd <- NULL
for (i in 1:5){
  i <- i/100
  cat(i)
  a <- sum(Cbd_FDR_MXRc$FDR_traditional_combined < i & Cbd_FDR_MXRc$Truth == T)
  b <- sum(Cbd_FDR_MXRc$Truth == T)
  
  sens_Cbd[i*100] <- (b-a)/b
  names(sens_Cbd)[i*100] <- i*100
}

spe <- rbind(spec, spec_Cbd)
sen <- rbind(sens, sens_Cbd)

sen[2,] - sen[1,]


############# PLOTS ################
#### Figure 1e,f ###


pdf(pdfname <- "Figure1.pdf", width=12, height=4)

layout(matrix(ncol=3, nrow=1, c(1:3)))
plot(y=1:dim(trad_FDR_MXR)[1],x=sort(trad_FDR_MXR$FDR_traditional), type="l", 
     ylab="recalled inter-links", xlab="target-decoy FDR", col="mediumseagreen", lwd=3, xlim=c(0,0.05), frame=F)
lines(y=1:dim(Cbd_FDR_MXRc)[1],x=sort(Cbd_FDR_MXRc$FDR_traditional_combined), type="l", col="orange", lwd=3)

legend("topleft", c("combined inter+intra FDR", "inter FDR"), col=c("orange", "mediumseagreen"), bty="n", lwd=3)


barplot(1-round(spe,3), ylim=c(0,0.3), ylab = "fraction false positive inter-links", xlab="target-decoy FDR-cutoff [%]", beside = T, main="false positives", col=rep(c("mediumseagreen", "orange"),5))
legend("topleft", c("combined inter+intra FDR", "inter FDR"), col=c("orange", "mediumseagreen"), bty="n", pch=15)

barplot(round(sen,3), ylim=c(0,0.55), ylab = "fraction false negative inter-links", xlab="target-decoy FDR-cutoff [%]", beside = T, main="false negatives",col=rep(c("mediumseagreen", "orange"),5))
legend("topleft", c("combined inter+intra FDR", "inter FDR"), col=c("orange", "mediumseagreen"), bty="n", pch=15)

graphics.off()
system(paste("open", pdfname))

#### Figure 1 bc ####

pdf(pdfname <- "Figure1_2.pdf", width=5, height=8)

layout(matrix(ncol=1,nrow=2, c(1,2)))
plot(density(XLs_ori_P$max_score[XLs_ori_P$Truth == "DECOY"]), col="tomato", lwd=3, xlim=c(0,30), xlab="-log10 score", main="target-decoy")
lines(density(XLs_ori_P$max_score[XLs_ori_P$Truth != "DECOY" & XLs_ori_P$Truth != "INTRALINK"]), col="steelblue2", lwd=3)
lines(density(XLs_ori_P$max_score[XLs_ori_P$Truth == F]), col="orange", lwd=3, xlim=c(0,30),ylim=c(0,.2), xlab="-log10 score", main="ground truth")
lines(density(XLs_ori_P$max_score[XLs_ori_P$Truth == T]), col="steelblue2", lwd=3)

graphics.off()
system(paste("open", pdfname))

