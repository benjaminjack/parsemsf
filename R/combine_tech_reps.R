#' Average top three peptides from multiple replicates
#'
#' @param df 
#'
#' @return
#' @export
#'
#' @examples
merge_top_peptides <- function(df, num_reps) {
  
  df %>%
    group_by(tech_rep, Sequence) %>%
    distinct(Area) %>% # Remove peptides with identical areas
    summarize(Area = sum(Area)) %>% # Sum across peptides with different charges
    top_n(3, desc(Area)) %>% # Find top 3 peptides
    ungroup() %>%
    group_by(Sequence) %>%
    mutate(n = n()) %>% # Count number of matched peptides
    ungroup() %>%
    filter(n >= num_reps) %>% # Remove any unmatched peptides
    summarize(area_mean = mean(Area), # Compute mean areas from top peptides
              area_sd = sd(Area),
              matched_peps = n()) -> matched_areas
  
  return(matched_areas)
}

#' Combine technical replicates using top 3 methods
#'
#' @param rep1 
#' @param rep2 
#'
#' @return
#' @export
#'
#' @examples
combine_tech_reps <- function(reps) {
  
  # A list to hold dataframes
  reps_df <- list()
  
  # Add an ID number for technical replicates
  for (i in 1:length(reps)) {
    reps_df[[i]] <- make_auc_table(reps[[i]]) %>% mutate(tech_rep = i)
  }
  
  combined <- bind_rows(reps_df)
  
  # How many technical replicates do we have?
  num_reps <- length(reps_df)
  
  combined %>%
    group_by(Proteins) %>%
    do(merge_top_peptides(., num_reps)) -> combined
  
  return(combined)
  
}