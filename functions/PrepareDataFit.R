
PrepareDataFit <- function(XLs_ori, SampSize){

  max_score <- -log10(apply(cbind(XLs_ori$n_score_a_MS2_MS3, XLs_ori$n_score_b_MS2_MS3),1, max))
  df_group <- cbind.data.frame(XLs_ori$LinkPos1,XLs_ori$protein_a_cln,  ## make a dataframe with the useful data
                               XLs_ori$LinkPos2,XLs_ori$protein_b_cln, max_score) ##
  
  names(df_group) <- c("K_alpha", "uniprot_alpha", "K_beta", "uniprot_beta", "max_score")
  
  AlphaProteins <- df_group$uniprot_alpha
  BetaProteins <-  df_group$uniprot_beta
  XLs <- cbind.data.frame(AlphaProteins, BetaProteins, max_score, XLs_ori$target_decoy, df_group$K_alpha, df_group$K_beta)
  
  names(XLs)[3] <- "Score"
  names(XLs)[5] <- "K_alpha"
  names(XLs)[6] <- "K_beta"
  
  XLs1 <- XLs[XLs_ori$protein_a_cln != XLs_ori$protein_b_cln,] ## Interlinks
  XLs2 <- XLs[XLs_ori$protein_a_cln == XLs_ori$protein_b_cln,] ## Intralinks

  ### estimated false positives ### take decoy distribution ###
  Ad <- sort(table(c(XLs$AlphaProteins[grepl(XLs$AlphaProteins, pattern = "#RND")], 
                     XLs$BetaProteins[grepl(XLs$BetaProteins, pattern = "#RND")])), decreasing=T)
  lam <- fitdistr(Ad, "Poisson")
  vecs <- rnorm(n = SampSize, mean=lam$estimate, sd = sqrt(lam$estimate))
  vecs[vecs<0] = 0
  vecs <- vecs/sum(vecs)
  
  ### simulated true positives ### take target entries
  XL_r <- XLs1[XLs1$`XLs_ori$target_decoy` == "target" & XLs1$Score>0,]
  A <- sort(table(c(XL_r$AlphaProteins, XL_r$BetaProteins)), decreasing=T)
  
  fr <- as.vector(A)
  p <- fr/sum(fr)
  ll <- function(s) sum(fr*(s*log(1:length(A))+log(sum(1/(1:length(A))^s))))
  lzipf <- function(s,N) -s*log(1:N)-log(sum(1/(1:N)^s))
  fit <- mle(ll,start=list(s=1))
  s.ll <- coef(fit)
  p2 <- exp(lzipf(s.ll, SampSize))
  
  return(list(p2,vecs))
}