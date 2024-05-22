

#### Figure 2 ####

setwd('N:/Interlink_subgroup/Scripts_functions')
sapply(paste('./Functions',list.files(pattern="[.]R$", path="./Functions"),sep='/'), source);

setwd("N:/Interlink_subgroup/Max_plates")

Truth <- read.csv(file = "all_plates_rerun_selectedfiles_Claasen_etAl_XlinkX.csv")

Truth <- Truth[!(grepl(Truth$protein_a, pattern="contam") | grepl(Truth$protein_b, pattern="contam")),] #& !(Truth$group_a != "" | Truth$group_b != "")),]

Truth <- Truth[!Truth$Truth == "CONTAMINANT",]

XLs_ori_1 <- Truth 
max_score <- -log10(apply(cbind(XLs_ori_1$n_score_a_MS2_MS3, XLs_ori_1$n_score_b_MS2_MS3),1, max))
XLs_ori_1 <- cbind.data.frame(XLs_ori_1, max_score)

XLs_ori_P <- XLs_ori_1
XLs_ori_P$ID <- rownames(XLs_ori_P)

Truth$ID <- rownames(Truth)

####### INTER-DEPENDENT #######

### fused search: groups
Links2_out <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_P, FUSED=T, Inter2=T, Grouped=T)
Links1_out <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_P, FUSED=T, Inter2=F, Grouped=T)
Old_MXR <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_P, FUSED=T, Inter2=F, Grouped=F)

### fused search: subgroup FDR
FDR_LI1_MXR <- MakeFDR2(sel=Links1_out, FDR_cut = 0.01)
FDR_LI2_MXR <- MakeFDR2(sel=Links2_out, FDR_cut = 0.01)
Old_FDR_MXR <- MakeFDR2(sel=Old_MXR, FDR_cut = 0.01)

### concatenated search: groups 
Links2_outc <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_P, FUSED=F, Inter2=T, Grouped=T)
Links1_outc <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_P, FUSED=F, Inter2=F, Grouped=T)
Old_MXRc <- SelectInterLinkDependentGroups(XLs_ori = XLs_ori_P, FUSED=F, Inter2=F, Grouped=F)

### concatenated search: subgroup FDR
FDR_LI1_MXRc <- MakeFDR2(sel=Links1_outc, FDR_cut = 0.01)
FDR_LI2_MXRc <- MakeFDR2(sel=Links2_outc, FDR_cut = 0.01)
Old_FDR_MXRc <- MakeFDR2(sel=Old_MXRc, FDR_cut = 0.01)


####### INTRA-DEPENDENT #######

### fused search: groups
MaxBI <- SelectSubset(XLs_ori = XLs_ori_P, BothIntras = T, Grouped = T, Concatenated = F)
MaxOI <- SelectSubset(XLs_ori = XLs_ori_P, BothIntras = F, Grouped = T, Concatenated = F)
MaxOld <- SelectSubset(XLs_ori = XLs_ori_P, BothIntras = F, Grouped = F, Concatenated = F)

### concatenated search: groups
MaxBIc <- SelectSubset(XLs_ori = XLs_ori_P, BothIntras = T, Grouped = T, Concatenated = T)
MaxOIc <- SelectSubset(XLs_ori = XLs_ori_P, BothIntras = F, Grouped = T, Concatenated = T)
MaxOldc <- SelectSubset(XLs_ori = XLs_ori_P, BothIntras = F, Grouped = F, Concatenated = T)

### fused search: subgroup FDR
FDR_MxR_both <- MakeFDR2(sel=MaxBI, FDR_cut = 0.01)
FDR_MxR_none <- MakeFDR2(sel=MaxOI, FDR_cut = 0.01)
FDR_MxR_old <- MakeFDR2(sel=MaxOld, FDR_cut = 0.01)

### concatenated search: subgroup FDR
FDR_MxR_bothc <- MakeFDR2(sel=MaxBIc, FDR_cut = 0.01)
FDR_MxR_nonec <- MakeFDR2(sel=MaxOIc, FDR_cut = 0.01)
FDR_MxR_oldc<- MakeFDR2(sel=MaxOldc, FDR_cut = 0.01)


####### HARMONIZE SUBGROUP FDRs #######

### INTER-DEPENDENT: fused 
INT_FDR_MXR <- HarmonizeSubgroupFDR(SG1=FDR_MxR_both, SG2 = FDR_MxR_none, FDR_append = "intra-dependent", XLs = XLs_ori_P, SG = T)
SG_trad_fus  <- HarmonizeSubgroupFDR(SG1=FDR_MxR_old, FDR_append = "traditional", XLs = XLs_ori_P, SG = F)

### INTRA-DEPENDENT: cocnatenated
INT_FDR_MXRc <- HarmonizeSubgroupFDR(SG1=FDR_MxR_bothc, SG2 = FDR_MxR_nonec, FDR_append = "intra-dependent-conc", XLs = XLs_ori_P, SG = T)
SG_trad_conc      <- HarmonizeSubgroupFDR(SG1=FDR_MxR_old, FDR_append = "traditional-conc", XLs = XLs_ori_P, SG = F)

### INTER-DEPENDENT: fused search 
Het_FDR_MXR  <-HarmonizeSubgroupFDR(SG1=FDR_LI1_MXR, SG2 = FDR_LI2_MXR, FDR_append = "inter-dependent", SG=T, XLs = XLs_ori_P)
trad_FDR_MXR <-HarmonizeSubgroupFDR(SG1=Old_FDR_MXR, FDR_append = "traditional", SG=F, XLs = XLs_ori_P)

### INTRA-DEPENDENT: concatenated
Het_FDR_MXRc  <-HarmonizeSubgroupFDR(SG1=FDR_LI1_MXRc, SG2 = FDR_LI2_MXRc, FDR_append = "inter-dependent_conc", SG=T, XLs = XLs_ori_P)
trad_FDR_MXRc <-HarmonizeSubgroupFDR(SG1=Old_FDR_MXRc, FDR_append = "traditional_conc", SG=F, XLs = XLs_ori_P)

####### MAKE GLOBAL FDRs from subgroup FDR #######

### inter-dependent # fused
ID_hSG <- cbind.data.frame(MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_MXR, removeDecoys = T, NO_SUBGROUPS = F)[[3]], 
                           MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_MXR, removeDecoys = T, NO_SUBGROUPS = F)[[4]])
names(ID_hSG) <- c("FDR_inter.dependent", "ID")

### inter-dependent # concatenated 
ID_hSGc <- cbind.data.frame(MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_MXRc, removeDecoys = T, NO_SUBGROUPS = F)[[3]], 
                            MakeFDR3(FDR_cut = 0.01, sel = Het_FDR_MXRc, removeDecoys = T, NO_SUBGROUPS = F)[[4]])
names(ID_hSGc) <- c("FDR_inter.dependent_conc", "ID")


### traditional 
ID_NG <- cbind.data.frame(MakeFDR3(FDR_cut = 0.01,  sel = trad_FDR_MXR,removeDecoys = T, NO_SUBGROUPS = F)[[3]], 
                          MakeFDR3(FDR_cut = 0.01,  sel = trad_FDR_MXR, removeDecoys = T, NO_SUBGROUPS = F)[[4]])
names(ID_NG) <- c("FDR_traditional", "ID")


### intra-dependent
ID_iSG <- cbind.data.frame(MakeFDR3(FDR_cut = 0.01, sel = INT_FDR_MXR, removeDecoys = T, NO_SUBGROUPS = F)[[3]], 
                           MakeFDR3(FDR_cut = 0.01, sel = INT_FDR_MXR, removeDecoys = T, NO_SUBGROUPS = F)[[4]])
names(ID_iSG) <- c("FDR_intra.dependent", "ID")

ID_iSGc <- cbind.data.frame(MakeFDR3(FDR_cut = 0.01, sel = INT_FDR_MXRc, removeDecoys = T, NO_SUBGROUPS = F)[[3]], 
                            MakeFDR3(FDR_cut = 0.01, sel = INT_FDR_MXRc, removeDecoys = T, NO_SUBGROUPS = F)[[4]])
names(ID_iSGc) <- c("FDR_intra.dependent_conc", "ID")

### merge IDs
id_list <- list(ID_hSG, ID_NG, ID_iSG, ID_iSGc, ID_hSGc)
Truth <- Reduce(function(x, y) merge(x, y, by = "ID"), id_list, init = Truth)



#### Function to calculate specificity (fraction false positives) and sensitivity (fraction false negatives) ####

calculate_spec_sens <- function(fdr_method, truth, cutoffs) {
  results <- list()
  for (cutoff in cutoffs) {
    a <- sum(Truth[[fdr_method]] < cutoff)
    b <- sum(Truth[[fdr_method]] < cutoff & Truth$Truth == FALSE)
    spec <- 1 - (b / a)
    sens <- (sum(Truth[[fdr_method]] < cutoff & Truth$Truth == TRUE)) / sum(Truth$Truth == "TRUE")
    results[[as.character(cutoff * 100)]] <- list(Specificity = spec, Sensitivity = sens)
  }
  return(results)
}

#### end function 



# Define cutoffs
cutoffs <- seq(0.005, 0.05, by = 0.005)

# Calculate Specificity and Sensitivity for each FDR method
spec_sens_results <- list(
  inter_intra_separate = calculate_spec_sens("FDR_traditional", Truth, cutoffs),
  Inter_dependent_fus = calculate_spec_sens("FDR_inter.dependent", Truth, cutoffs),
  Inter_dependent_conc = calculate_spec_sens("FDR_inter.dependent_conc", Truth, cutoffs),
  Intra_dependent_fus = calculate_spec_sens("FDR_intra.dependent", Truth, cutoffs),
  Intra_dependent_conc = calculate_spec_sens("FDR_intra.dependent_conc", Truth, cutoffs)
)
spec_results <- lapply(spec_sens_results, function(x) sapply(x, function(y) y$Specificity))
sens_results <- lapply(spec_sens_results, function(x) sapply(x, function(y) y$Sensitivity))

spec_df <- as.data.frame(spec_results)
sens_df <- as.data.frame(sens_results)

fpr_results <- 1 - spec_df
fnr_results <- 1 - sens_df

###### Plotting Figure 2b,c and 4 a,b #####
par(mfrow=c(1, 2))  

barplot(t(fpr_results), beside=TRUE, col=c('darkgrey','#fdae61','#d7191c','#abd9e9','#2c7bb6') , ylim=c(0, 0.18), 
        main="false positives", xlab="target-decoy FDR cut-off [%]", ylab="fraction false positives",
        legend.text = colnames(fpr_results), args.legend = list(x="topleft", bty="n"))
grid(NULL,NULL)
barplot(t(fnr_results), beside=TRUE, col=c('darkgrey','#fdae61','#d7191c','#abd9e9','#2c7bb6'), ylim=c(0, 0.65), 
        main="false negatives", xlab="target-decoy FDR cut-off [%]", ylab="fraction false negatives",
        legend.text = colnames(fnr_results), args.legend = list(x="topright", bty="n"))

