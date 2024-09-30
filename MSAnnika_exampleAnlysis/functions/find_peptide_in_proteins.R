### this function finds a peptide and lysine position in a protein

find_peptide_in_proteins <- function(peptide, sequences) {

  peptide_clean <- str_remove_all(peptide, "\\[|\\]")
  
  k_pos_peptide <- stri_locate_first_fixed(peptide, "[K]")
  k_pos_peptide <- k_pos_peptide[1]
  
  results <- logical(length(sequences))
  k_positions_in_protein <- integer(length(sequences))
  
  for (i in seq_along(sequences)) {
    seq_string <- sequences[i]
    match_pos <- stri_locate_first_fixed(seq_string, peptide_clean)
    
    if (!is.na(match_pos)) {
      results[i] <- TRUE
      k_positions_in_protein[i] <- match_pos[1] + k_pos_peptide - 1
    }
  }
  return(list(
    found = results,
    k_positions = k_positions_in_protein[k_positions_in_protein >0]
  ))
}
