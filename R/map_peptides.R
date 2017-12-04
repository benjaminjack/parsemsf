
#' Map peptides to their locations within a protein
#'
#' Takes a ThermoFisher MSF file and finds the location of each peptide within its corresponding protein sequence. In cases where a single peptide maps to multiple locations within a protein sequence, only the first location is reported. If a peptide maps ambiguously to multiple proteins, all locations are reported with data from each peptide-protein combination on a separate row.
#'
#' @inheritParams make_pep_table
#'
#' @return A dataframe containing start and stop positions (relative to the parent protein sequence) for each peptide in the database.
#'
#' \item{peptide_id}{a unique peptide ID}
#' \item{spectrum_id}{a unique spectrum ID}
#' \item{protein_id}{unique protein group ID to which this peptide maps}
#' \item{protein_desc}{protein description from reference database used to assign peptides to protein groups, parsed according to \code{prot_regex}}
#' \item{peptide_sequence}{amino acid sequence (does not show post-translational modifications)}
#' \item{pep_score}{PEP score}
#' \item{q_value}{Q-value score}
#' \item{protein_sequence}{parent protein sequence}
#' \item{start}{start position of peptide within protein sequence}
#' \item{end}{end position of peptide within protein sequence}
#'
#' @export
#'
#' @examples
#' map_peptides(parsemsf_example("test_db.msf"))
map_peptides <- function(msf_file, min_conf = "High", prot_regex = "") {

  pep_table <- make_pep_table(msf_file, min_conf, collapse = FALSE) %>%
    rename_(.dots = setNames(list(~sequence), c("peptide_sequence")))

  my_db <- DBI::dbConnect(RSQLite::SQLite(), msf_file)

  prots <- tbl(my_db, "Proteins") %>%
    select_(~ProteinID, ~Sequence) %>%
    collect()

  pep_table <- inner_join(pep_table, prots, by = c("protein_id" = "ProteinID")) %>%
    rename_(.dots = setNames(list(~Sequence), c("protein_sequence"))) %>%
    rowwise() %>%
    mutate_(.dots = setNames(list(~str_locate(protein_sequence, peptide_sequence)[1],
                                  ~str_locate(protein_sequence, peptide_sequence)[2]),
                             c("start",
                               "end")))
  DBI::dbDisconnect(my_db)

  return(pep_table)
}
