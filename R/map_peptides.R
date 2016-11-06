
#' Map peptides to their locations within a protein
#'
#' Takes a Thermo MSF file and finds the location of each peptide within its corresponding protein sequence. In cases where a single peptide maps to multiple locations within a protein sequence, only the first location is reported.
#'
#' @inheritParams make_pep_table
#'
#' @return A dataframe containing start and stop positions (relative to the parent protein sequence) for each peptide in the database.
#'
#' \item{PeptideID}{unique peptide ID}
#' \item{SpectrumID}{unique spectrum ID}
#' \item{Pep_seq}{peptide sequence}
#' \item{Proteins}{protein description}
#' \item{ProteinID}{unique proteinID}
#' \item{PEP}{?}
#' \item{q-Value}{?}
#' \item{Prot_seq}{parent protein sequence}
#' \item{start}{start position of peptide within protein sequence}
#' \item{end}{end position of peptide within protein sequence}
#'
#' @export
#'
#' @examples
#' map_peptides(parsemsf_example("test_db.msf"))
map_peptides <- function(msf_file, min_conf = "High", prot_regex = "") {

  pep_table <- make_pep_table(msf_file, min_conf, collapse = FALSE) %>%
    rename_(.dots = setNames(list(~Sequence), c("Pep_seq")))

  my_db <- src_sqlite(msf_file)

  prots <- tbl(my_db, "Proteins") %>%
    select_(~ProteinID, ~Sequence) %>%
    collect()

  pep_table <- inner_join(pep_table, prots, by = "ProteinID") %>%
    rename_(.dots = setNames(list(~Sequence), c("Prot_seq"))) %>%
    rowwise() %>%
    mutate_(.dots = setNames(list(~str_locate(Prot_seq, Pep_seq)[1],
                                  ~str_locate(Prot_seq, Pep_seq)[2]),
                             c("start",
                               "end")))

  return(pep_table)
}
