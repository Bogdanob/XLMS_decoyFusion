
SelectInterLinkDependentGroups2 <- function(XLs_ori, FUSED=T, Inter2=T, Grouped=T){
  
 # XLs_ori <- XLs; FUSED = T; Grouped=T; Inter2=T
  max_score <- XLs_ori$Best.CSM.Score
  df_group <- cbind.data.frame(XLs_ori$In.protein.A,XLs_ori$Accession.A,  ## make a dataframe with the useful data
                               XLs_ori$In.protein.B,XLs_ori$Accession.B, max_score) ##  
  names(df_group) <- c("K_alpha", "uniprot_alpha", "K_beta", "uniprot_beta", "max_score")
  
  AlphaProteins <- df_group$uniprot_alpha
  BetaProteins <-  df_group$uniprot_beta
  XLs <- cbind.data.frame(AlphaProteins, BetaProteins, max_score, XLs_ori$targettarget,XLs_ori$targetdecoy, XLs_ori$decoydecoy, df_group$K_alpha, df_group$K_beta, XLs_ori$ID)
  
  names(XLs)[3] <- "max_score"
  names(XLs)[4] <- "targettarget"
  names(XLs)[5] <- "targetdecoy"
  names(XLs)[6] <- "decoydecoy"
  names(XLs)[7] <- "K_alpha"
  names(XLs)[8] <- "K_beta"
  names(XLs)[9] <- "ID"
  ### Start making  LysLys counts per PPI
  
  Tabulated<- table(apply(cbind.data.frame(XLs$AlphaProteins,XLs$BetaProteins), 1, function(x) paste(sort(x), collapse="_"))) ## how many Lys-Lys per PPI
  Ints1  <- apply(cbind.data.frame(XLs$AlphaProteins,XLs$BetaProteins), 1, function(x) paste(sort(x), collapse="_")) ## unique PPI identifier
  
  I1 <- apply(cbind.data.frame(XLs$AlphaProteins,XLs$K_alpha),1, function(x) paste(x ,collapse="_K_"))
  I2 <- apply(cbind.data.frame(XLs$BetaProteins,XLs$K_beta),1, function(x) paste(x ,collapse="_K_"))
  Ints2  <- apply(cbind.data.frame(I1,I2), 1, function(x) paste(sort(x), collapse="_")) ## unique PPI identifier
  Ints2 <- str_remove_all(Ints2, " ")
  LysLysPPI <- aggregate(Ints2, by=list(Ints1), unique)
  LysLysPPI$x <- lapply(LysLysPPI$x, paste, collapse=" ")
  
  # Extract cross-linked lysines as a nested list
  KKs <- lapply(LysLysPPI$x, str_extract_all, "K_[0-9]+")

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
  
  if (Grouped ==T & Inter2 ==F){
    Output  <- XLs[Ints1_Bin1,]
  }
  
  if (Grouped ==T & Inter2 ==T){
    Output  <- XLs[Ints1_Bin2,]
  }
  
  if (Grouped ==F){
    Output  <- XLs
  }
  
  return(Output)
}
