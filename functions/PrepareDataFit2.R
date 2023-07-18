library(MASS)
library(stringr)

PrepareDataFit2 <- function(XLs_ori, SampSize){
  
   SampSize <- 4860
  
   XLs_ori <- XLs_ori_CW
  max_score <- -log10(apply(cbind(XLs_ori$n_score_a_MS2_MS3, XLs_ori$n_score_b_MS2_MS3),1, max))
  df_group <- cbind.data.frame(XLs_ori$LinkPos1,XLs_ori$gene_a,  ## make a dataframe with the useful data
                               XLs_ori$LinkPos2,XLs_ori$gene_b, max_score) ##
  
  names(df_group) <- c("K_alpha", "uniprot_alpha", "K_beta", "uniprot_beta", "max_score")
  
  AlphaProteins <- df_group$uniprot_alpha
  BetaProteins <-  df_group$uniprot_beta
  XLs <- cbind.data.frame(AlphaProteins, BetaProteins, -log10(max_score), XLs_ori$target_decoy, df_group$K_alpha, df_group$K_beta)
  
  names(XLs)[3] <- "Score"
  names(XLs)[5] <- "K_alpha"
  names(XLs)[6] <- "K_beta"
  
  XLs1 <- XLs[XLs_ori$protein_a_cln != XLs_ori$protein_b_cln,] ## Interlinks

  Adb4  <- sort(table(apply(cbind(XLs_ori$gene_a[XLs$`XLs_ori$target_decoy`=="decoy"],XLs_ori$gene_b[XLs$`XLs_ori$target_decoy`=="decoy"]), 1, function(x) paste(sort(x), collapse=","))), decreasing=T)
  
  
  
  
  CollapseFPs <- function(i){
    selvec <- names(sort(table(unlist(lapply(str_split(string=names(Adb4),pattern=","),"[[",1))), decreasing=T)[i])
    Ad <- sort(as.vector(Adb4[(grepl(pattern = selvec,names(Adb4)))]), decreasing=T)  
    return(Ad)
    }
  
  C1 <- c(CollapseFPs(i))
  C3<-NULL
  for (i in 1:100) {
    C2 <- C1
    C1 <- (c(CollapseFPs(i)))
    C3 <- c(C3,C1)
    }
  
  
  #Ad <- sort(table(c(XL_d$AlphaProteins[grepl(XL_d$AlphaProteins, pattern = "#RND")], 
   #                  XL_d$BetaProteins[grepl(XL_d$BetaProteins, pattern = "#RND")])), decreasing=T)
  lam <- fitdistr(sort(C3, decreasing=T), "Poisson")
  vecs <- rnorm(n = SampSize, mean=lam$estimate, sd = sqrt(lam$estimate))
  vecs[vecs<0] = 0
  vecs <- vecs/sum(vecs)
  


  ### estimated true positives ### targets ### 
  

  Ab4 <- table(apply(cbind(XLs1$AlphaProteins[XLs1$`XLs_ori$target_decoy`=="target"],XLs1$BetaProteins[XLs1$`XLs_ori$target_decoy`=="target"]), 1, function(x) paste(sort(x), collapse=",")) )

  sel <- names(sort(table(unlist(lapply(str_split(string=names(Ab4),pattern=","),"[[",1))), decreasing=T))
  
 
  
  lengis <- NULL
  crossis <- NULL
  for (i in 1:length(sel)){
    A <- sort(as.vector(Ab4[(grepl(pattern = sel[i],names(Ab4)))]), decreasing=T)
    lengis[i] <- length(A)
    crossis[i] <- sum(A)
    }
  
  lm3 <- lm(log2(lengis)~ log2(crossis))
  
  
  A <- sort(as.vector(Ab4[(grepl(pattern = sel[1],names(Ab4)))]), decreasing=T)
  
  fr <- as.vector(A)
  p <- fr/sum(fr)
  ll <- function(s) sum(fr*(s*log(1:length(A))+log(sum(1/(1:length(A))^s))))
  fit <- mle(ll,start=list(s=1))
  s.ll <- coef(fit)
  p2 <- exp(lzipf(s.ll, SampSize))
  
  return(list(p2,vecs, lm3, s.ll))
}