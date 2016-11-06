
#' Make a data frame of peptides from a Thermo MSF file
#'
#' @description Extracts amino acid sequences (without post-translational modifications), protein matches, and quality scores.
#' @param msf_file A file path to a thermo MSF file.
#' @param min_conf "High", "Medium", or "Low". The mininum peptide confidence level to retrieve from MSF file.
#' @param prot_regex Regular expression where the first group matches a protein name or ID from the protein description. Regex must contain ONE group. The protein description is typically generated from a fasta reference file that was used for the database search.
#' @param collapse If TRUE, peptides that match to multiple protein sequences are collapsed into a single row with multiple protein descriptions and protein IDs in the \code{Proteins} and \code{ProteinID} columns separated by semi-colons (";").
#'
#' @return A data frame of all peptides above the confidence cut-off from a thermo MSF file.
#'
#' \item{PeptideID}{a unique peptide ID}
#' \item{SpectrumID}{a unique spectrum ID}
#' \item{Sequence}{amino acid sequence (does not show post-translational modifications)}
#' \item{Proteins}{protein name from protein description in reference database, parsed according to \code{prot_regex}}
#' \item{ProteinID}{unique protein ID to which this peptide maps}
#' \item{PEP}{?}
#' \item{q-Value}{?}
#'
#' @export
#'
#' @examples
#' # Read from a path
#'
#' make_pep_table(parsemsf_example("test_db.msf"))
#'
#' # Retrieve full protein description
#'
#' make_pep_table(parsemsf_example("test_db.msf"), prot_regex = "")
#'
#' # ...which is also equivalent to...
#'
#' make_pep_table(parsemsf_example("test_db.msf"), prot_regex = "^(.+)$")
#'
make_pep_table <- function(msf_file,
                           min_conf = "High",
                           prot_regex = "^>([a-zA-Z0-9._]+)\\b",
                           collapse = TRUE) {

  confidence <- switch(min_conf,
                       High = 3,
                       Medium = 2,
                       Low = 1,
                       3)

  if (prot_regex == "") {
    prot_regex = "^(.+)$"
  }

  # Access MSF database file
  my_db <- src_sqlite(msf_file)

  # This table maps PeptideIDs to ProteinIDs
  PeptidesProteins <- tbl(my_db, "PeptidesProteins") %>%
    select_(~PeptideID, ~ProteinID)

  # Here are the actual peptide sequences with corresponding SpectrumIDs
  Peptides <- tbl(my_db, "Peptides") %>%
    select_(~PeptideID, ~SpectrumID, ~ConfidenceLevel, ~Sequence) %>%
    filter_(~ ConfidenceLevel >= confidence)

  # Protein descriptions are contained in this table
  ProteinAnnotations <- tbl(my_db, "ProteinAnnotations") %>%
    select_(~ProteinID, ~Description)

  # Lazy-eval formula for protein name matching
  prot_match <- lazyeval::interp(~str_match(Description, prot_regex)[,2], prot_regex = prot_regex)

  # Build a peptide table
  pep_table <- inner_join(PeptidesProteins, Peptides, by = c("PeptideID" = "PeptideID")) %>%
    inner_join(ProteinAnnotations, by = "ProteinID") %>%
    collect(n = Inf) %>%
    # Extract protein ID; assumes that protein ID is immediately after ">" and ends with a space
    mutate_(.dots = setNames(list(prot_match), c("Proteins"))) %>%
    group_by_(~PeptideID)
  # Collapse peptides that map to multiple proteins into a single row
  if (collapse == TRUE) {
    pep_table <- pep_table %>%
      summarize_(.dots = setNames(list(~unique(SpectrumID),
                                       ~unique(Sequence),
                                       ~paste(Proteins, sep = "; "),
                                       ~paste(ProteinID, sep = "; ")),
                                  c("SpectrumID",
                                    "Sequence",
                                    "Proteins",
                                    "ProteinID")))
  }
  pep_table <- pep_table %>% select_(~PeptideID, ~SpectrumID, ~Sequence, ~Proteins, ~ProteinID)
  # Append custom fields
  CustomFields <- tbl(my_db, "CustomDataFields")

  # Grab custom peptide data using SQL because of automatic sqlite typing
  # This assumes that FieldValue will always be a number
  CustomPeptides <- tbl(my_db, sql("SELECT FieldID, PeptideID, CAST(FieldValue as REAL) AS FieldValue FROM CustomDataPeptides"))

  # Spread custom fields as separate columns
  custom_data <- left_join(CustomPeptides, CustomFields, by = "FieldID") %>%
    select_(~PeptideID, ~DisplayName, ~FieldValue) %>%
    collect(n = Inf) %>%
    spread_("DisplayName", "FieldValue")

  # Join to peptide table
  pep_table <- left_join(pep_table, custom_data, by = "PeptideID")

  return(pep_table)

}
