
### intra-dependent grouping depending on whether concatenated or fused database is used

MakeSubgroupFDR <- function(XLs_ori, Concatenated = F){
  
  if (Concatenated == T){
    df_group1 <-  subset(XLs_ori$protein_a_cln != XLs_ori$protein_b_cln, x=XLs_ori)## Interlinks
    df_group2 <-  subset(XLs_ori$protein_a_cln == XLs_ori$protein_b_cln, x=XLs_ori)## Intralinks
    
    HaveIntra <- unique(c(df_group2$protein_a_cln, df_group2$protein_b_cln)) ### proteins with an intra-link
    NoIntra_sep <- df_group1[(!df_group1$protein_a_cln %in% HaveIntra) & (!df_group1$protein_b_cln %in% HaveIntra),] ### dataframe of links, proteins have no intra-link
    BothIntra_sep <- df_group1[(df_group1$protein_a_cln %in% HaveIntra) & (df_group1$protein_b_cln %in% HaveIntra),] ### dataframe of links, both proteins have intra-link
    AtLeastOneIntra_sep <- df_group1[(df_group1$protein_a_cln %in% HaveIntra) | (df_group1$protein_b_cln %in% HaveIntra),] ### dataframe of links, at least one of the proteins have intra-link
    OneIntra_sep <- df_group1[((df_group1$protein_a_cln %in% HaveIntra) | (df_group1$protein_b_cln %in% HaveIntra)) &! ### dataframe of links, one protein has intra-link
                                ((df_group1$protein_a_cln %in% HaveIntra) & (df_group1$protein_b_cln %in% HaveIntra)),] ### dataframe of links, <2 proteins have intra-link
    OneIntraNoIntra_sep <- rbind(OneIntra_sep, NoIntra_sep) ### dataframe of links, <2 proteins have intra-link
    
    BI_sep <- table(BothIntra_sep$target_decoy == "decoy")
    NI_sep <- table(NoIntra_sep$target_decoy == "decoy")
    OI_sep <- table(OneIntra_sep$target_decoy == "decoy")
    OINI_sep <- table(OneIntraNoIntra_sep$target_decoy == "decoy")
    
    if (is.na(BI_sep[2])){
      BI_sep[2] <- 0
    }
    return(matrix(c(BI_sep, OINI_sep), ncol=2,nrow=2))  
  }
  
  if (Concatenated == F){   #### target decoy fusion #### 
    df_group1 <-  subset(XLs_ori$protein_a_cln != XLs_ori$protein_b_cln, x=XLs_ori)## Interlinks
    df_group2 <-  subset(XLs_ori$protein_a_cln == XLs_ori$protein_b_cln, x=XLs_ori)## Intralinks
    
    HaveIntra <- unique(c(df_group2$gene_a, df_group2$gene_b)) ### proteins with an intra-link
    NoIntra_fus <- df_group1[(!df_group1$gene_a %in% HaveIntra) & (!df_group1$gene_b %in% HaveIntra),]
    BothIntra_fus <- df_group1[(df_group1$gene_a %in% HaveIntra) & (df_group1$gene_b %in% HaveIntra),]
    AtLeastOneIntra <- df_group1[(df_group1$gene_a %in% HaveIntra) | (df_group1$gene_b %in% HaveIntra),]
    OneIntra_fus <- df_group1[((df_group1$gene_a %in% HaveIntra) | (df_group1$gene_b %in% HaveIntra)) &! 
                                ((df_group1$gene_a %in% HaveIntra) & (df_group1$gene_b %in% HaveIntra)),]
    OneIntraNoIntra_fus <- rbind(OneIntra_fus, NoIntra_fus)
    
    
    BI_fus <- table(BothIntra_fus$target_decoy == "decoy")
    NI_fus <- table(NoIntra_fus$target_decoy == "decoy")
    OI_fus <- table(OneIntra_fus$target_decoy == "decoy")
    OINI_fus <- table(OneIntraNoIntra_fus$target_decoy == "decoy")
    
    return(matrix(c(BI_fus, OINI_fus), ncol=2,nrow=2))  
  }
}
