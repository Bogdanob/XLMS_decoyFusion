
SelectSubset <- function(XLs_ori = XLs_ori, BothIntras = F, Grouped = F, Concatenated = F){
  
#  XLs_ori <- XLs_ori_P
#  Grouped=T
#  Concatenated=F
# BothIntras = F
  
  if (Grouped ==T & Concatenated ==F){
    df_group1 <-  subset(XLs_ori$gene_a != XLs_ori$gene_b, x=XLs_ori)## Interlinks
    df_group2 <-  subset(XLs_ori$gene_b == XLs_ori$gene_a, x=XLs_ori)## Intralinks
    
    HaveIntra <- unique(c(df_group2$gene_a, df_group2$gene_b)) ### proteins with an intra-link
    NoIntra_fus <- df_group1[(!df_group1$gene_a %in% HaveIntra) & (!df_group1$gene_b %in% HaveIntra),]
    BothIntra_fus <- df_group1[(df_group1$gene_a %in% HaveIntra) & (df_group1$gene_b %in% HaveIntra),]
    AtLeastOneIntra <- df_group1[(df_group1$gene_a %in% HaveIntra) | (df_group1$gene_b %in% HaveIntra),]
    OneIntra_fus <- df_group1[((df_group1$gene_a %in% HaveIntra) | (df_group1$gene_b %in% HaveIntra)) &! 
                                ((df_group1$gene_a %in% HaveIntra) & (df_group1$gene_b %in% HaveIntra)),]
    OneIntraNoIntra_fus <- rbind(OneIntra_fus, NoIntra_fus)
    
    if (BothIntras == T){
      Output <- BothIntra_fus
    }  
    
    if (BothIntras == F){
      Output <- OneIntraNoIntra_fus
    } 
  }  
  
  if (Grouped ==T & Concatenated ==T){
    df_group1 <-  subset(XLs_ori$protein_a_cln != XLs_ori$protein_b_cln, x=XLs_ori)## Interlinks
    df_group2 <-  subset(XLs_ori$protein_a_cln == XLs_ori$protein_b_cln, x=XLs_ori)## Intralinks
    
    HaveIntra <- unique(c(df_group2$protein_a_cln, df_group2$protein_b_cln)) ### proteins with an intra-link
    NoIntra_fus <- df_group1[(!df_group1$protein_a_cln %in% HaveIntra) & (!df_group1$protein_b_cln %in% HaveIntra),]
    BothIntra_fus <- df_group1[(df_group1$protein_a_cln %in% HaveIntra) & (df_group1$protein_b_cln %in% HaveIntra),]
    AtLeastOneIntra <- df_group1[(df_group1$protein_a_cln %in% HaveIntra) | (df_group1$protein_b_cln %in% HaveIntra),]
    OneIntra_fus <- df_group1[((df_group1$protein_a_cln %in% HaveIntra) | (df_group1$protein_b_cln %in% HaveIntra)) &! 
                                ((df_group1$protein_a_cln %in% HaveIntra) & (df_group1$protein_b_cln %in% HaveIntra)),]
    OneIntraNoIntra_fus <- rbind(OneIntra_fus, NoIntra_fus)
    
    if (BothIntras == T){
      Output <- BothIntra_fus
    }  
    
    if (BothIntras == F){
      Output <- OneIntraNoIntra_fus
    } 
  }
  
  
  
  
  if (Grouped == F){   #### target decoy concatenated #### 
    df_group1 <-  subset(XLs_ori$protein_a_cln != XLs_ori$protein_b_cln, x=XLs_ori)## Interlinks
    df_group2 <-  subset(XLs_ori$protein_a_cln == XLs_ori$protein_b_cln, x=XLs_ori)## Intralinks
    
    Output <- df_group1
    
    
  }
  return(Output)
}
