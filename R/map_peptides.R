
#' Map peptides to their locations within a protein
#'
#' Takes a thermo MSF file and uses `stringr` to find the location of peptides within their corresponding protein sequences.
#'
#' @param msf_file
#' @param min_conf
#'
#' @return A dataframe containing start and stop positions (relative to the parent protein sequence) for each peptide in the database.
#' @export
#'
#' @examples
#' map_peptides("mythermofile.msf")
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
