
make_pep_table <- function(msf_file, min_conf = "High") {
  
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
  
  # Append custom fields
  CustomFields <- tbl(my_db, "CustomDataFields")
  
  # Grab custom peptide data using SQL because of automatic sqlite typing
  # This assumes that FieldValue will always be a number
  CustomPeptides <- tbl(my_db, sql("SELECT FieldID, PeptideID, CAST(FieldValue as REAL) AS FieldValue FROM CustomDataPeptides"))
  
  # Spread custom fields as separate columns
  custom_data <- left_join(CustomPeptides, CustomFields) %>%
    select(PeptideID, DisplayName, FieldValue) %>%
    collect() %>%
    spread(DisplayName, FieldValue)
  
  # Join to peptide table
  pep_table <- left_join(pep_table, custom_data)
  
  return(pep_table)

}