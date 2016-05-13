#' Make a data frame of peptide scores
#'
#' Makes a data frame containing scores (PEP, Q-value) for each peptides in a thermo MSF file.
#'
#' @param msf_file A file path to a thermo MSF file.
#' @param min_conf "High", "Medium", or "Low". The mininum peptide confidence level to retrieve from MSF file.
#'
#' @return A data frame containing PEP and Q-value for peptides in an MSF file.
#' @export
#'
#' @examples
#' make_score_table("mythermofile.msf")
make_score_table <- function(msf_file, min_conf = "High") {

  # Access MSF database file
  my_db <- src_sqlite(msf_file)

  PeptideScores <- tbl(my_db, "PeptideScores")

  ProcessingNodeScores <- tbl(my_db, "ProcessingNodeScores")

  score_table <- left_join(PeptideScores, ProcessingNodeScores, by="ScoreID") %>%
    select(PeptideID, ScoreName, ScoreValue) %>%
    collect() %>%
    spread(ScoreName, ScoreValue)

  return(score_table)

}
