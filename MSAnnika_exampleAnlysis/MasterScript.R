##################
#### MSAnnika #### 
##################

# Exemplary analysis of context-sensitive FDR for data searched by MSAnnika #
### Steps ###
## 1) read CSM data from MSAnnika
## 2) perform CSM FDR on MSAnnika data
## 3) read MSAnnika cross-link level data and subset to those XLs surviving CSM FDR
## 4) match decoy sequences to proteins (not outputted by MSAnnika)...takes a while
## 5) inter-dependent subgrouping XL-level
## 6) perform FDR using PEPs and integrate into final list
## 7) generate some plots highlighting the benefit


#### Load necessary library ####
library(Biostrings)
library(stringr)
library(seqinr)
library(stringi)
library(data.table)


########## 1 ##########
#### Load CSM data
CSMs <- read.csv("N:/Interlink_subgroup/revision files/Cong_HEK_DSBSO/DSBSO_HEK_MSAnnika_CSMs.txt", sep="\t")

#### restrict analysis to inter-links
CSMs <- subset(CSMs,CSMs$Crosslink.Type == "Inter")

### define targets and decoys
CSMs$targettarget <- CSMs$Accession.A != "" & CSMs$Accession.B != ""
CSMs$decoydecoy <- CSMs$Accession.A == "" & CSMs$Accession.B == ""
CSMs$targetdecoy <- (CSMs$Accession.A == "" | CSMs$Accession.B == "") & !CSMs$decoydecoy


### give a unique ID to CSMs
CSMs$ID <- rownames(CSMs)

########## 2 ##########
#### Make CSM FDR, when NO_subgroups == T, makes standard FDR, FDR_cut as set is not important but returns in the first list entry the score at which cut-off should be performed to limit decoys to set percentage
CSMFDR<- MakePEP_FDR3(sel=CSMs, MSAnnika = T, FDR_cut = 0.2, NO_SUBGROUPS = T)

# Filter CSM data by Interlink FDR, arbitrarily set to 10 %. 
CSMs_filt <- subset(CSMFDR[[5]], CSMFDR[[5]]$FDR < 0.1)

## get XL identifier
CSM_sequence <- paste(CSMs_filt$Sequence.A,CSMs_filt$Crosslinker.Position.A, 
                      CSMs_filt$Sequence.B,CSMs_filt$Crosslinker.Position.B, sep="-")

########## 3 ##########
# now read the XL-file
XLs <- read.csv("N:/Interlink_subgroup/revision files/Cong_HEK_DSBSO/DSBSO_HEK_MSAnnika_Crosslinks.txt", sep="\t")

# get XL-identifier 
XL_posA <- unlist(lapply(str_locate_all(pattern = "\\[", XLs$Sequence.A), "[[", 1))
XL_posB <- unlist(lapply(str_locate_all(pattern = "\\[", XLs$Sequence.B), "[[", 1))

XL_sequence <- str_remove_all(pattern = "\\[", string=paste(XLs$Sequence.A,XL_posA,  XLs$Sequence.B,XL_posB, sep="-"))
XL_sequence <- str_remove_all(pattern = "\\]", string=XL_sequence)

## subset XL-list by surviving CSMs

XLs<- subset(XLs, XL_sequence %in% CSM_sequence)
fasta <- read.fasta("N:/Interlink_subgroup/revision files/Cong_HEK_DSBSO/20220323_CW_DSBSO_Paper_HEK_db.fasta", as.string = T, seqtype = "AA")

########## 4 ##########
# make decoy fasta

reversed_fasta_dt <- data.table(
  name = names(fasta),
  sequence = sapply(fasta, function(entry) reverse_sequence(entry[[1]]))
)
# Extract sequences
sequences <- reversed_fasta_dt$sequence

# remember decoys
XLs$targettarget <- XLs$Accession.A != "" & XLs$Accession.B != ""
XLs$decoydecoy <- XLs$Accession.A == "" & XLs$Accession.B == ""
XLs$targetdecoy <- (XLs$Accession.A == "" | XLs$Accession.B == "")& (!XLs$decoydecoy)

### LOOKUP decoy matching IDs
# Preprocess peptides, restrict decoys
preprocessed_peptides_A <- XLs$Sequence.A[which(XLs$Accession.A == "")]
preprocessed_peptides_B <- XLs$Sequence.B[which(XLs$Accession.B == "")]

# Process peptides (match to decoy IDs #)
results_A <- process_peptides(preprocessed_peptides_A, sequences, reversed_fasta_dt$name)
results_B <- process_peptides(preprocessed_peptides_B, sequences, reversed_fasta_dt$name)

decoy_indices_A <- which(XLs$Accession.A == "")
decoy_indices_B <- which(XLs$Accession.B == "")

## write the decoy protein into XL-file
XLs$Accession.A[decoy_indices_A] <- unlist(lapply(str_extract_all(pattern = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", string = results_A$names_list), paste0,collapse=";"))
XLs$Accession.B[decoy_indices_B] <- unlist(lapply(str_extract_all(pattern = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", string = results_B$names_list), paste0,collapse=";"))

XLs$In.protein.A[decoy_indices_A] <- results_A$k_list
XLs$In.protein.B[decoy_indices_B] <- results_B$k_list

XLs <- subset(XLs, !XLs$Accession.A == XLs$Accession.B)
XLs$ID <-  row.names(XLs)

########## 6 ##########
Links2_out_CW <- SelectInterLinkDependentGroups2(XLs_ori = XLs, FUSED=T, Inter2=T, Grouped=T) ## context-rich
Links1_out_CW <- SelectInterLinkDependentGroups2(XLs_ori = XLs, FUSED=T, Inter2=F, Grouped=T) ## context-poor
Old_CW        <- SelectInterLinkDependentGroups2(XLs_ori = XLs, FUSED=T, Inter2=F, Grouped=F) ## standard approach

PEP_Links1_out <- MakePEP3(sel = Links1_out_CW, MSAnnika = T) ## PEPs context-poor
PEP_Links2_out <- MakePEP3(sel = Links2_out_CW, MSAnnika = T) ## PEPs context-rich
PEP_OLD_out <- MakePEP3(sel = Old_CW, MSAnnika = T) ## standard approach

## now make the FDRs based on the PEPs
Interdependent_Groups <-MakePEP_FDR3(FDR_cut = 0.01, sel = rbind(PEP_Links1_out, PEP_Links2_out), 
                                     removeDecoys = T, MSAnnika = T)[[5]]
No_Groups<- MakePEP_FDR3(FDR_cut = 0.01, sel = PEP_OLD_out, 
                      removeDecoys = T, MSAnnika = T, NO_SUBGROUPS = T)[[5]]


########## 7 ##########
#### PLOTS ####
pdf(pdfname <- "N:/Interlink_subgroup/ReviewerFigureMSAnnika.pdf", width=18, height=6)
layout(matrix(ncol=3,nrow=1, c(1:3)))

plot(sort(CSMFDR[[5]]$FDR[CSMFDR[[5]]$targettarget == T], decreasing = F), 1:length(CSMFDR[[5]]$FDR[CSMFDR[[5]]$targettarget == T]), type="l", 
     xlim=c(0,0.2), ylim=c(0,21000), xlab="target-decoy FDR [%]", ylab="number of inter-CSMs")
abline(v=0.1, col="darkgrey", lty=3)
grid(NULL,NULL)


barplot(matrix(ncol=2, nrow=2, c(table(Links1_out_CW$targettarget),table(Links2_out_CW$targettarget))), 
        col=c("red", "lightgrey"), names.arg = c("context-poor", "context-rich"), ylab="inter-link count")

plot(y=1:length(No_Groups$FDR),sort(No_Groups$FDR), type="l", xlim=c(0,0.13), ylim=c(0,9000), ylab="number of inter-links (Residue pairs)", xlab="target-decoy FDR cut-off [%]")
lines(y=1:length(Interdependent_Groups$FDR),sort(Interdependent_Groups$FDR), col="red")

#sum(Interdependent_Groups$FDR < 0.01)/sum(No_Groups$FDR < 0.01)
abline(v=0.01, lty=3)

graphics.off()
system(paste("open", pdfname))
