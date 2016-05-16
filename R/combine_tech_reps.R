#' Average areas of top three peptides from multiple technical replicates
#'
#' Takes a data frame containing peptide information for a single protein group and several technical replicates and averages the most abundant peptides (by area) across all replicates. By default, this function takes the averages across the top three peptides. This function will also optionally match peptides across replicates so only matching peptides are averaged together.
#'
#' @param df A data frame generated from \code{\link{make_auc_table}} for a single protein group with the column \code{tech_rep} to indicate the technical replicate
#' @param num_reps numeric Number of technical replicates being combined
#' @param match_peps boolean Should we only quantitate based on the top three peptides present in all replicates?
#'
#' @return A data frame with the columns \code{area_mean}, \code{area_sd}, \code{peps_per_tech_rep}, which corresponds to the average area, the standard deviation of the areas, and the number of peptides that were averaged together divided by the number of replicates. This dataframe corresponds to a single protein group.
#' @export
#'
#' @examples
#' merge_top_peptides(df, 2, match_peps = T)
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

#' Combine technical replicates and quantitate proteins
#'
#' Takes a list of thermo MSF files, parses and combines them into a single data frame, and computes areas for each protein group based on the top 3 method of quantitation.
#'
#' @param reps Vector. List of thermo MSF file names
#' @param normalize Boolean. Should we normalize peptide areas for technical replicate to the total areas in a given replicate?
#' @param match_peps Boolean. Should we quantitate only on matching peptides across technical replicates?
#' @param relabel Named vector for relabeling protein groups. Names correspond to a pattern or string to match (i.e. the name or ID of a protein group), and values correspond to the new value (i.e. new protein group name).
#'
#' @return A data frame with the columns \code{area_mean}, \code{area_sd}, \code{peps_per_tech_rep}, which corresponds to the average area, the standard deviation of the areas, and the number of peptides that were averaged together divided by the number of replicates for each protein group. This dataframe contains all protein groups.
#' @export
#'
#' @examples
#' combine_tech_reps(c("rep1.msf", "rep2.msf"), relabel = c("NP_12345.1" = "NP_1000.1"))
combine_tech_reps <- function(reps, normalize = TRUE, match_peps = TRUE, relabel = c()) {

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
  if (length(relabel) > 0) {
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
