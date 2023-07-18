GeneNameExtractXlinkX <- function(VirTable){
  
  # VirTable <- Data
  # VirTable$gene_a[!grepl(x=VirTable$protein_a, pattern="lcl")] <- str_extract(string = VirTable$protein_a[!grepl(x=VirTable$protein_a, pattern="lcl")], pattern="GN=[^ ]{2,9}")
  #  VirTable$gene_b[!grepl(x=VirTable$protein_b, pattern="lcl")] <- str_extract(string = VirTable$protein_b[!grepl(x=VirTable$protein_b, pattern="lcl")], pattern="GN=[^ ]{2,9}")
  
  temp <- VirTable
  
  
  VirTable$gene_a[grepl(x=VirTable$protein_a, pattern="lcl")] <- str_extract(string = VirTable$protein_a[grepl(x=VirTable$protein_a, pattern="lcl")], pattern="protein=[^ ]{2,9}")
  VirTable$gene_b[grepl(x=VirTable$protein_b, pattern="lcl")] <- str_extract(string = VirTable$protein_b[grepl(x=VirTable$protein_b, pattern="lcl")], pattern="protein=[^ ]{2,9}")
  
  VirTable$gene_a[grepl(x=temp$protein_a, pattern="lcl")] <- str_extract(string=str_extract(string = VirTable$gene_a, pattern="=[A-Za-z0-9.]{2,9}"), pattern= "[A-Za-z0-9.]{2,9}")[grepl(x=temp$protein_a, pattern="lcl")]   
  VirTable$gene_b[grepl(x=temp$protein_b, pattern="lcl")] <- str_extract(string=str_extract(string = VirTable$gene_b, pattern="=[A-Za-z0-9.]{2,9}"), pattern= "[A-Za-z0-9.]{2,9}")[grepl(x=temp$protein_b, pattern="lcl")]  
  return(VirTable)  
}
