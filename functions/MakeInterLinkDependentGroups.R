
MakeInterLinkDependentGroups <- function(XLs_ori=XLs_ori, FUSED=T){
  max_score <- apply(cbind(XLs_ori$n_score_a_MS2_MS3, XLs_ori$n_score_b_MS2_MS3),1, max)
  
  
  #### change here for target decoy (fusion/concatenated)
  if (FUSED==F){
    df_group <- cbind.data.frame(XLs_ori$LinkPos1,XLs_ori$protein_a_cln,  ## make a dataframe with the useful data
                                 XLs_ori$LinkPos2,XLs_ori$protein_b_cln, max_score) ##  
  }
  if (FUSED==T){
    df_group <- cbind.data.frame(XLs_ori$LinkPos1,XLs_ori$gene_a,  ## make a dataframe with the useful data
                                 XLs_ori$LinkPos2,XLs_ori$gene_b, max_score) ##  
  }
  
  names(df_group) <- c("K_alpha", "uniprot_alpha", "K_beta", "uniprot_beta", "max_score")
  
  #MaxResID <- aggregate(x = df_group$score, by=list(ResID), max) # take the maximum score from all CSMs for a specifc lys-lys ID
  
  AlphaProteins <- df_group$uniprot_alpha
  BetaProteins <-  df_group$uniprot_beta
  
  
  XLs <- cbind.data.frame(AlphaProteins, BetaProteins, -log10(max_score), XLs_ori$target_decoy, df_group$K_alpha, df_group$K_beta)
  
  names(XLs)[3] <- "Score"
  names(XLs)[5] <- "K_alpha"
  names(XLs)[6] <- "K_beta"
  
  XLs1 <- XLs[XLs_ori$protein_a_cln != XLs_ori$protein_b_cln,] ## Interlinks
  XLs2 <- XLs[XLs_ori$protein_a_cln == XLs_ori$protein_b_cln,] ## Intralinks
  
  ### Start trying to make Bins here
  
  Tabulated<- table(apply(cbind.data.frame(XLs1$AlphaProteins,XLs1$BetaProteins), 1, function(x) paste(sort(x), collapse="_"))) ## how many Lys-Lys per PPI
  Ints1  <- apply(cbind.data.frame(XLs1$AlphaProteins,XLs1$BetaProteins), 1, function(x) paste(sort(x), collapse="_")) ## unique PPI identifier
  
  I1 <- apply(cbind.data.frame(XLs1$AlphaProteins,XLs1$K_alpha),1, function(x) paste(x ,collapse="_K_"))
  I2 <- apply(cbind.data.frame(XLs1$BetaProteins,XLs1$K_beta),1, function(x) paste(x ,collapse="_K_"))
  
  Ints2  <- apply(cbind.data.frame(I1,I2), 1, function(x) paste(sort(x), collapse="_")) ## unique PPI identifier
  
  Ints2 <- str_remove_all(Ints2, " ")
  
  LysLysPPI <- aggregate(Ints2, by=list(Ints1), unique)
  
  LysLysPPI$x <- lapply(LysLysPPI$x, paste, collapse=" ")
  
  # Extract cross-linked lysines as a nested list
  KKs <- lapply(LysLysPPI$x, str_extract_all, "K_[0-9]*")
  
  # Flatten the list
  KKs_flat <- unlist(KKs, recursive = FALSE)
  
  # Create a list of even and odd lysines for each interaction
  KKs_even_odd <- lapply(KKs_flat, function(x) {
    KKs_even <- x[c(FALSE, TRUE)]
    KKs_odd <- x[c(TRUE, FALSE)]
    list(KKs_even = KKs_even, KKs_odd = KKs_odd)
  })
  
  # Check the output
  KKs_counts <- lapply(KKs_even_odd, function(x) {
    KKs_even_count <- table(unlist(x$KKs_even))
    KKs_odd_count <- table(unlist(x$KKs_odd))
    list(KKs_even_count = KKs_even_count, KKs_odd_count = KKs_odd_count)
  })
  
  KKs_bins <- lapply(KKs_counts, function(x) {
    ifelse(length(x$KKs_odd_count) > 1 & length(x$KKs_even_count) > 1, "Both > 1", 
           ifelse(length(x$KKs_odd_count) == 1 | length(x$KKs_even_count) == 1, "at least One = 1"))
  })
  
  MasterTable <- cbind.data.frame(unlist(KKs_bins), LysLysPPI$Group.1, unlist(LysLysPPI$x))
  
  Ints1_Bin1 <- Ints1 %in% MasterTable$`LysLysPPI$Group.1`[unlist(KKs_bins)=="at least One = 1"]
  Ints1_Bin2 <- Ints1 %in% MasterTable$`LysLysPPI$Group.1`[unlist(KKs_bins)=="Both > 1"]
  
  
  return(matrix(ncol=2, nrow=2, c(table(XLs1$`XLs_ori$target_decoy`[Ints1_Bin1] == "decoy"), table(XLs1$`XLs_ori$target_decoy`[Ints1_Bin2] == "decoy"))))
  #  return(matrix(ncol=2, nrow=1, table(XLs1$`XLs_ori$target_decoy`[Ints1_Bin2] == "decoy")
  
}