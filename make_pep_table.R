
make_pep_table <- function(msf_file, min_conf = 3) {
  
  confidence <- switch(min_conf,
                       High = 3,
                       Medium = 2,
                       Low = 1,
                       3)
  
  # Access MSF database file
  my_db <- src_sqlite(msf_file)
  
  # This table maps PeptideIDs to ProteinIDs
  PeptidesProteins <- tbl(my_db, "PeptidesProteins") %>%
    select(PeptideID, ProteinID)
  
  # Here are the actual peptide sequences with corresponding SpectrumIDs
  Peptides <- tbl(my_db, "Peptides") %>%
    select(PeptideID, SpectrumID, ConfidenceLevel, Sequence) %>%
    filter(ConfidenceLevel >= confidence)
  
  # Protein descriptions are contained in this table
  ProteinAnnotations <- tbl(my_db, "ProteinAnnotations") %>%
    select(ProteinID, Description)
  
  # Build a peptide table
  pep_table <- inner_join(PeptidesProteins, Peptides, by = c("PeptideID" = "PeptideID")) %>% 
    inner_join(ProteinAnnotations) %>%
    collect() %>%
    mutate(Proteins = str_match(Description, "^>([a-zA-Z0-9._]+)\\b")[,2]) %>%
    group_by(PeptideID) %>%
    summarize(SpectrumID = unique(SpectrumID), 
              Sequence = unique(Sequence), 
              Proteins = paste(Proteins, collapse = "; "))
  
  return(pep_table)

}