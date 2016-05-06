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
combine_tech_reps <- function(reps, normalize = TRUE, match_peps = TRUE) {
  
  # A list to hold dataframes
  reps_df <- list()
  
  # Add an ID number for technical replicates
  for (i in 1:length(reps)) {
    reps_df[[i]] <- make_auc_table(reps[[i]]) %>% 
      mutate(tech_rep = i)
  }
  
  
  # Rename some protein groups
  # This should probably not be hard-coded but I haven't come up with a smart way to do it automatically
  bind_rows(reps_df) %>%
    mutate(Proteins = str_replace(Proteins, 'NP_041997.1; NP_041998.1', 'NP_041998.1')) %>%
    mutate(Proteins = str_replace(Proteins, 'NP_041975.1; NP_041977.1', 'NP_041975.1')) %>%
    mutate(Proteins = str_replace(Proteins, 'NP_041997.1', 'NP_041998.1')) %>%
    mutate(Proteins = str_replace(Proteins, 'NP_041977.1', 'NP_041975.1')) -> combined
  
  # How many technical replicates do we have?
  num_reps <- length(reps_df)
  
  # Check if we should normalize to total area for a given replicate
  # This accounts for any variability in how the sample was injected
  if (normalize == TRUE) {
    combined %>%
      group_by(tech_rep) %>%
      mutate(total_area = sum(Area, na.rm = TRUE), Area = Area/total_area) %>%
      ungroup() %>%
      group_by(Proteins) %>%
      do(merge_top_peptides(., num_reps, match_peps = match_peps)) -> combined
  } else {
    combined %>%
      group_by(Proteins) %>%
      do(merge_top_peptides(., num_reps, match_peps = match_peps)) -> combined
  }
  
  return(combined)
  
}