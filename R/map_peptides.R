
#' Map peptides to their locations within a protein
#'
#' Takes a Thermo MSF file and finds the location of each peptide within its corresponding protein sequence. In cases where a single peptide maps to multiple locations within a protein sequence, only the first location is reported.
#'
#' @inheritParams make_pep_table
#'
#' @return A dataframe containing start and stop positions (relative to the parent protein sequence) for each peptide in the database.
#' @export
#'
#' @examples
#' # map_peptides("mythermofile.msf")
map_peptides <- function(msf_file, min_conf = "High", prot_regex = "") {

  pep_table <- make_pep_table(msf_file, min_conf, collapse = FALSE) %>%
    rename(Pep_seq = Sequence)

  my_db <- src_sqlite(msf_file)

  prots <- tbl(my_db, "Proteins") %>%
    select(ProteinID, Sequence) %>%
    collect()

  pep_table <- inner_join(pep_table, prots, by = "ProteinID") %>%
    rename(Prot_seq = Sequence) %>%
    rowwise() %>%
    mutate(start = str_locate(Prot_seq, Pep_seq)[1],
           end = str_locate(Prot_seq, Pep_seq)[2])
}
