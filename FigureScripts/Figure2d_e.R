
### Figure 2 d,e ###
### ROC curve ###
library(ROCit)

### Max ### ground-truth comparison
setwd('N:/Interlink_subgroup/Scripts_functions')
sapply(paste('./Functions',list.files(pattern="[.]R$", path="./Functions"),sep='/'), source);

setwd("N:/Interlink_subgroup/Max_plates")

Truth <- read.csv(file = "Truth__matched_contaminantsincluded.csv")
dim(Truth)
### how to deal with contaminants # remove
Truth <- Truth[!(grepl(Truth$protein_a, pattern="contam") | grepl(Truth$protein_b, pattern="contam")),] #& !(Truth$group_a != "" | Truth$group_b != "")),]




XLs_ori_1 <- read.csv("N:/Interlink_subgroup/Max_plates//SS_unique_lys_crosslink_plate1.csv") 
max_score <- -log10(apply(cbind(XLs_ori_1$n_score_a_MS2_MS3, XLs_ori_1$n_score_b_MS2_MS3),1, max))
XLs_ori_1 <- cbind.data.frame(XLs_ori_1, max_score)

XLs_ori_2 <- read.csv("N:/Interlink_subgroup/Max_plates//SS_unique_lys_crosslink_plate2.csv") 
max_score <- -log10(apply(cbind(XLs_ori_2$n_score_a_MS2_MS3, XLs_ori_2$n_score_b_MS2_MS3),1, max))
XLs_ori_2 <- cbind.data.frame(XLs_ori_2, max_score)

XLs_ori_3 <- read.csv("N:/Interlink_subgroup/Max_plates//SS_unique_lys_crosslink_plate3.csv") 
max_score <- -log10(apply(cbind(XLs_ori_3$n_score_a_MS2_MS3, XLs_ori_3$n_score_b_MS2_MS3),1, max))
XLs_ori_3 <- cbind.data.frame(XLs_ori_3, max_score)

XLs_ori_P <- rbind(XLs_ori_1, XLs_ori_2, XLs_ori_3)


XLs_ori_P$ID <- rownames(XLs_ori_P)

XLs_ori_P$gene_a[(is.na(XLs_ori_P$gene_a))] <- str_extract(XLs_ori_P$protein_a_cln[(is.na(XLs_ori_P$gene_a))], pattern = '.*NO ')
XLs_ori_P$gene_a[(is.na(XLs_ori_P$gene_a))] <- str_extract(XLs_ori_P$protein_a_cln[(is.na(XLs_ori_P$gene_a))], pattern = '.*NO$')
XLs_ori_P$gene_a <- str_remove_all(pattern = "###RND###", XLs_ori_P$gene_a)

XLs_ori_P$gene_b[(is.na(XLs_ori_P$gene_b))] <- str_extract(XLs_ori_P$protein_b_cln[(is.na(XLs_ori_P$gene_b))], pattern = '.*NO ')
XLs_ori_P$gene_b[(is.na(XLs_ori_P$gene_b))] <- str_extract(XLs_ori_P$protein_b_cln[(is.na(XLs_ori_P$gene_b))], pattern = '.*NO$')
XLs_ori_P$gene_b <- str_remove_all(pattern = "###RND###", XLs_ori_P$gene_b)



## deal with contaminants ## remove
Con <- (grepl(XLs_ori_P$protein_a, pattern="contam") & grepl(XLs_ori_P$protein_b, pattern="contam")) 
XLs_ori_P <- subset(XLs_ori_P, !Con) #| XLs_ori_P$target_decoy == "decoy")

### inter-dependent groups ### fused search 
Links2_out <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_P, FUSED=T, Inter2=T, Grouped=T)
Links1_out <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_P, FUSED=T, Inter2=F, Grouped=T)
Old_MXR <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_P, FUSED=T, Inter2=F, Grouped=F)

### make subgroup FDR ### inter-link dependent ### fused
FDR_LI1_MXR <- MakeFDR2(sel=Links1_out, FDR_cut = 0.01)
FDR_LI2_MXR <- MakeFDR2(sel=Links2_out, FDR_cut = 0.01)
Old_FDR_MXR <- MakeFDR2(sel=Old_MXR, FDR_cut = 0.01)


### inter-dependent groups ### concatenated search 
Links2_outc <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_P, FUSED=F, Inter2=T, Grouped=T)
Links1_outc <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_P, FUSED=F, Inter2=F, Grouped=T)
Old_MXRc <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_P, FUSED=F, Inter2=F, Grouped=F)

#Links2_outc <- subset(Links2_outc, !Con | Links2_outc$target_decoy == "decoy")
#Links1_outc <- subset(Links1_outc, !Con | Links1_outc$target_decoy == "decoy")
#Old_MXR <- subset(Old_MXR, !Con | Old_MXR$target_decoy == "decoy")

### make subgroup FDR ### inter-link dependent ### fused
FDR_LI1_MXRc <- MakeFDR2(sel=Links1_outc, FDR_cut = 0.01)
FDR_LI2_MXRc <- MakeFDR2(sel=Links2_outc, FDR_cut = 0.01)
Old_FDR_MXRc <- MakeFDR2(sel=Old_MXRc, FDR_cut = 0.01)


### intra-link dependent Groups ### fused search
MaxBI <- SelectSubset(XLs_ori = XLs_ori_P, BothIntras = T, Grouped = T, Concatenated = F)
MaxOI <- SelectSubset(XLs_ori = XLs_ori_P, BothIntras = F, Grouped = T, Concatenated = F)
MaxOld <- SelectSubset(XLs_ori = XLs_ori_P, BothIntras = F, Grouped = F, Concatenated = F)


### intra-link dependent groups ### concatenated search 
MaxBIc <- SelectSubset(XLs_ori = XLs_ori_P, BothIntras = T, Grouped = T, Concatenated = T)
MaxOIc <- SelectSubset(XLs_ori = XLs_ori_P, BothIntras = F, Grouped = T, Concatenated = T)
MaxOldc <- SelectSubset(XLs_ori = XLs_ori_P, BothIntras = F, Grouped = F, Concatenated = T)


### make subgroup FDR ### intra-dependent ### fused
FDR_MxR_both <- MakeFDR2(sel=MaxBI, FDR_cut = 0.01)
FDR_MxR_none <- MakeFDR2(sel=MaxOI, FDR_cut = 0.01)
FDR_MxR_old <- MakeFDR2(sel=MaxOld, FDR_cut = 0.01)

INT_FDR_MXR <- HarmonizeSubgroupFDR(SG1=FDR_MxR_both, SG2 = FDR_MxR_none, FDR_append = "intra-dependent", XLs = XLs_ori_P, SG = T)
SG_trad_fus      <- HarmonizeSubgroupFDR(SG1=FDR_MxR_old, FDR_append = "traditional", XLs = XLs_ori_P, SG = F)

### make subgroup FDR ### intra-dependent ### concatenated
FDR_MxR_bothc <- MakeFDR2(sel=MaxBIc, FDR_cut = 0.01)
FDR_MxR_nonec <- MakeFDR2(sel=MaxOIc, FDR_cut = 0.01)
FDR_MxR_oldc<- MakeFDR2(sel=MaxOldc, FDR_cut = 0.01)

INT_FDR_MXRc <- HarmonizeSubgroupFDR(SG1=FDR_MxR_bothc, SG2 = FDR_MxR_nonec, FDR_append = "intra-dependent-conc", XLs = XLs_ori_P, SG = T)
SG_trad_conc      <- HarmonizeSubgroupFDR(SG1=FDR_MxR_old, FDR_append = "traditional-conc", XLs = XLs_ori_P, SG = F)

### inter-link dependent groups ### fused search 
Het_FDR_MXR  <-HarmonizeSubgroupFDR(SG1=FDR_LI1_MXR, SG2 = FDR_LI2_MXR, FDR_append = "inter-dependent", SG=T, XLs = XLs_ori_P)
trad_FDR_MXR <-HarmonizeSubgroupFDR(SG1=Old_FDR_MXR, FDR_append = "traditional", SG=F, XLs = XLs_ori_P)

### inter-link dependent groups ### concatenated search 
Het_FDR_MXRc  <-HarmonizeSubgroupFDR(SG1=FDR_LI1_MXRc, SG2 = FDR_LI2_MXRc, FDR_append = "inter-dependent_conc", SG=T, XLs = XLs_ori_P)
trad_FDR_MXRc <-HarmonizeSubgroupFDR(SG1=Old_FDR_MXRc, FDR_append = "traditional_conc", SG=F, XLs = XLs_ori_P)


Truth$ID %in% MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_MXR, removeDecoys = T)[[4]]
ID_hSG <- cbind.data.frame(MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_MXR, removeDecoys = T)[[3]], MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_MXR, removeDecoys = T)[[4]])
names(ID_hSG) <- c("FDR_inter.dependent", "ID")
ID_hSGc <- cbind.data.frame(MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_MXRc, removeDecoys = T)[[3]], MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_MXRc, removeDecoys = T)[[4]])
names(ID_hSGc) <- c("FDR_inter.dependent_conc", "ID")
ID_NG <- cbind.data.frame(MakeFDR3(FDR_cut = 0.01,  sel = trad_FDR_MXR,removeDecoys = T)[[3]], MakeFDR3(FDR_cut = 0.01, sel = trad_FDR_MXR, removeDecoys = T)[[4]])
names(ID_NG) <- c("FDR_traditional", "ID")
ID_iSG <- cbind.data.frame(MakeFDR3(FDR_cut = 0.01, sel = INT_FDR_MXR, removeDecoys = T)[[3]], MakeFDR3(FDR_cut = 0.01, sel = INT_FDR_MXR, removeDecoys = T)[[4]])
names(ID_iSG) <- c("FDR_intra.dependent", "ID")
ID_iSGc <- cbind.data.frame(MakeFDR3(FDR_cut = 0.01, sel = INT_FDR_MXRc, removeDecoys = T)[[3]], MakeFDR3(FDR_cut = 0.01, sel = INT_FDR_MXRc, removeDecoys = T)[[4]])
names(ID_iSGc) <- c("FDR_intra.dependent_conc", "ID")

View(Truth)

Truth2 <- merge(Truth, ID_hSG, by="ID")
Truth3 <- merge(Truth2, ID_NG, by="ID")
Truth4 <- merge(Truth3, ID_iSG, by="ID")
Truth5 <- merge(Truth4, ID_iSGc, by="ID")
Truth6 <- merge(Truth5, ID_hSGc, by="ID")
Truth<-Truth6


CompareFDR <- function(cutoff){
# cutoff = 0.007
  
  ToMap_trad <- as.matrix(rev(table(Truth$Truth[Truth$FDR_traditional < cutoff])))
  ToMap_inter <- as.matrix(rev(table(Truth$Truth[Truth$FDR_inter.dependent < cutoff])))
  ToMap_intra <- as.matrix(rev(table(Truth$Truth[Truth$FDR_intra.dependent < cutoff])))
  ToMap_intra_conc <- as.matrix(rev(table(Truth$Truth[Truth$FDR_intra.dependent_conc < cutoff])))
  ToMap_inter_conc <- as.matrix(rev(table(Truth$Truth[Truth$FDR_inter.dependent_conc < cutoff])))
  
  return(as.matrix(cbind(ToMap_trad, ToMap_inter, ToMap_intra, ToMap_intra_conc, ToMap_inter_conc)))
}

data005 <- CompareFDR(0.006)



empFDR_trad<-NULL
empFDR_hSG<-NULL
empFDR_iSG<-NULL
empFDR_iSGc<-NULL
empFDR_hSGc<-NULL
FDR_tar<-NULL
for (i in c(6:200)){
cat(i)
  
   data <- CompareFDR(i/1000)
  empFDR_trad[i] <- data[2,1]/colSums(data)[1]
  empFDR_hSG[i] <- data[2,2]/colSums(data)[2]
  empFDR_iSG[i] <- data[2,3]/colSums(data)[3]
  empFDR_iSGc[i] <- data[2,4]/colSums(data)[4]
  empFDR_hSGc[i] <- data[2,5]/colSums(data)[5]
  
  FDR_tar[i] <-  i/1000

  }
pdf(pdfname <- "Figure2de.pdf", width=15, height=5)

layout(matrix(ncol=3, nrow=1, c(1:3)))

plot(FDR_tar, empFDR_trad, type="l", col="darkgrey", lwd=3, xlab="target FDR", ylab="empirical FDR", xlim=c(0,0.16), ylim=c(0,0.16))
lines(FDR_tar, empFDR_hSG, type="l", col="mediumseagreen", lwd=3)
lines(FDR_tar, empFDR_iSG, type="l", col="orange", lwd=3)
lines(FDR_tar, empFDR_iSGc, type="l", col="orange", lwd=3, lty=3)
lines(FDR_tar, empFDR_hSGc, type="l", col="mediumseagreen", lwd=3, lty=3)



abline(a=0,b=1, lty=2, lwd=3, col="lightgrey")
#segments(x0 = 0.01,x1=0.01, y0=0, y1=0.01, lty=2, lwd=3)
#segments(x0 = 0.0,x1=0.01, y0=0.01, y1=0.01, lty=2, lwd=3)
grid(NULL,NULL)

legend("bottomright", col = c("darkgrey","mediumseagreen","orange", "mediumseagreen", "orange"),title = "FDR classifier",
       c("no subgroups", "inter-dependent (fused)",
         "intra-dependent (fused)", "inter-dependent (concatenated)",
         "intra-dependent (concatenated)"), lwd = c(3,3,3,3,3),lty=c(1,1,1,3,3), cex=1)


rocit_trad <- rocit(score = 1/Truth$FDR_traditional, 
                    class = Truth$Truth, negref = F ,
                    method = "emp")
rocit_Inter.dep <- rocit(score = 1/Truth$FDR_inter.dependent, 
                         class = Truth$Truth,negref=F, 
                         method = "emp")
rocit_intra.dep <- rocit(score = 1/Truth$FDR_intra.dependent, 
                         class = Truth$Truth,negref = F ,
                         method = "emp")
ChanceClassifier <- rnorm(n = length(Truth$FDR_traditional))
rocit_chance <- rocit(score = ChanceClassifier, 
                      class = Truth$Truth,negref = F ,
                      method = "emp")
rocit_intra.dep_conc <- rocit(score = 1/Truth$FDR_intra.dependent_conc, class = Truth$Truth, negref=F, method="emp")
rocit_inter.dep_conc <- rocit(score = 1/Truth$FDR_inter.dependent_conc, class = Truth$Truth, negref=F, method="emp")






plot(rocit_trad, col = c("darkgrey"), 
     legend = FALSE, YIndex = FALSE)
lines(rocit_Inter.dep$TPR~rocit_Inter.dep$FPR, 
      col = "mediumseagreen", lwd = 2)
lines(rocit_intra.dep$TPR~rocit_intra.dep$FPR, 
      col = "orange", lwd = 2)
lines(rocit_chance$TPR~rocit_chance$FPR, 
      col = 6, lwd = 2)
lines(rocit_intra.dep_conc$TPR~rocit_intra.dep_conc$FPR, 
      col = "orange", lwd = 2, lty=3)
lines(rocit_inter.dep_conc$TPR~rocit_inter.dep_conc$FPR, 
      col = "mediumseagreen", lwd = 2, lty=3)

legend("bottomright", col = c("darkgrey","mediumseagreen","orange", "mediumseagreen", "orange"),title = "FDR classifier",
       c("no subgroups", "inter-dependent (fused)",
         "intra-dependent (fused)", "inter-dependent (concatenated)",
         "intra-dependent (concatenated)"), lwd = c(3,3,3,3,3),lty=c(1,1,1,3,3), cex=1)


grid(NULL,NULL)

graphics.off()
system(paste("open",pdfname))







