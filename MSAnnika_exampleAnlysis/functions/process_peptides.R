process_peptides <- function(peptides, sequences, names) {
  # peptides <- preprocessed_peptides_A[1] 
  #  names <- reversed_fasta_dt$name
  names_list <- c(NULL)
  k_list <- c(NULL)
  
  for (i in seq_along(peptides)) {
    peptide <- peptides[i]
    # Perform the search
    results <- find_peptide_in_proteins(peptide, sequences)
    results_names <- results$found
    # Collect matching identifiers
    matching_ids <- names[results_names]
    names_list[i] <- paste0(matching_ids, collapse = ";")
    k_list[i] <- paste0(results$k_positions, collapse = ";")
    cat(i / length(peptides) * 100, "% complete\r")
  }
  
  return(cbind.data.frame(names_list,k_list))
}