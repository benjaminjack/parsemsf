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