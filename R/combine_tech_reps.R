#' Average top three peptides from multiple replicates
#'
#' @param df 
#'
#' @return
#' @export
#'
#' @examples
merge_top_peptides <- function(df, num_reps, match_peps = TRUE) {
  
  df %>%
    group_by(tech_rep, Sequence) %>%
    distinct(Area) %>% # Remove peptides with identical areas
    summarize(Area = sum(Area, na.rm = T)) %>% # Sum across peptides with different charges
    group_by(tech_rep) %>%
    top_n(3, Area) %>% # Find top 3 peptides
    ungroup() -> df
  
  if (match_peps == TRUE) {
    df %>% group_by(Sequence) %>%
      mutate(n = n()) %>% # Count number of matched peptides
      ungroup() %>%
      filter(n >= num_reps) -> df # Remove any unmatched peptides
  }
  
  df %>%
    summarize(area_mean = mean(Area, na.rm = T), # Compute mean areas from top peptides
              area_sd = sd(Area, na.rm = T),
              peps_per_tech_rep = n()/num_reps) -> matched_areas
  
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
combine_tech_reps <- function(reps, normalize = TRUE, match_peps = TRUE, relabel = FALSE) {
  
  # A list to hold dataframes
  reps_df <- list()
  
  # How many technical replicates do we have?
  num_reps <- length(reps)
  
  # Add an ID number for technical replicates
  for (i in 1:length(reps)) {
    message(paste("Now processing: ", reps[[i]]))
    reps_df[[i]] <- make_auc_table(reps[[i]]) %>% 
      mutate(tech_rep = i)
  }
  
  # Combine into single dataframe with all technical replicates
  combined <- bind_rows(reps_df)
  
  # Rename some protein groups
  if (relabel != FALSE) {
    combined %<>% mutate(Proteins = str_replace_all(Proteins, relabel))
  }
  
  # Check if we should normalize to total area for a given replicate
  # This accounts for any variability in how the sample was injected
  message("Quantitating...")
  if (normalize == TRUE) {
    combined %<>%
      group_by(tech_rep) %>%
      mutate(total_area = sum(Area, na.rm = TRUE), Area = Area/total_area) %>%
      ungroup()
  } 
  
  # Quantitate using top three most abundant areas
  combined %<>%
    group_by(Proteins) %>%
    do(merge_top_peptides(., num_reps, match_peps = match_peps))
  
  return(combined)
  
}