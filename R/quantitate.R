#' Average areas of top three peptides from multiple technical replicates
#'
#' Takes a data frame containing peptide information for a single protein group and several technical replicates and averages the most abundant peptides (by area) across all replicates. By default, this function takes the averages across the top three peptides. This function will also optionally match peptides across replicates so only matching peptides are averaged together.
#'
#' @param df A data frame generated from \code{\link{make_area_table}} for a single protein group with the column \code{tech_rep} to indicate the technical replicate
#' @param num_reps numeric Number of technical replicates being combined
#' @param match_peps boolean Should we only quantitate based on peptides present in all replicates?
#'
#' @return A data frame that corresponds to a single protein group.
#'
#' \item{area_mean}{average peptide area}
#' \item{area_sd}{peptide area standard deviation}
#' \item{peps_per_rep}{number of peptides per technical replicate used to calculate area_mean and area_sd. This is typically 3, but may be less.}
#' @export
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{merge_top_peptides(df, 2, match_peps = T)}
merge_top_peptides <- function(df, num_reps, match_peps = TRUE) {

  df %>%
    group_by_(~tech_rep, ~sequence) %>%
    distinct_(.dots = list(~area)) %>% # Remove peptides with identical areas
    summarize_(.dots = setNames(list(~sum(area, na.rm = T)),
                                c("area"))) %>% # Sum across peptides with different charges
    group_by_(~tech_rep) %>%
    filter_(~ min_rank(desc(area)) <= 3) %>% # Find top 3 peptides
    ungroup() -> df

  if (match_peps == TRUE) {
    df %>% group_by_(~sequence) %>%
      mutate_(.dots = setNames(list(~n()), c("n"))) %>% # Count number of matched peptides
      ungroup() %>%
      filter_(~ n >= num_reps) -> df # Remove any unmatched peptides
  }

  dots <- setNames(list(~mean(area, na.rm = T),
                        ~sd(area, na.rm = T),
                        ~ n()/num_reps),
                   c("area_mean",
                     "area_sd",
                     "peps_per_rep"))

  df %>%
    summarize_(.dots = dots) -> matched_areas # Compute mean areas from top peptides

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
#' @return A data frame containing area information for all proteins.
#'
#' \item{protein_desc}{protein description}
#' \item{area_mean}{average peptide area}
#' \item{area_sd}{peptide area standard deviation}
#' \item{peps_per_rep}{Number of peptides per technical replicate used to calculate \code{area_mean} and \code{area_sd}. This is typically 3 peptides, but may be less.}
#'
#' @export
#'
#' @examples
#' quantitate(c(parsemsf_example("test_db.msf"),
#'              parsemsf_example("test_db2.msf")),
#'            relabel = c("NP_12345.1" = "NP_1000.1"))
quantitate <- function(reps, normalize = T, match_peps = T, relabel = c()) {

  # A list to hold dataframes
  reps_df <- list()

  # How many technical replicates do we have?
  num_reps <- length(reps)

  # Add an ID number for technical replicates
  for (i in 1:length(reps)) {
    message(paste("Now processing: ", reps[[i]]))
    reps_df[[i]] <- make_area_table(reps[[i]]) %>%
      mutate_(.dots = setNames(list(i), c("tech_rep")))
  }

  # Combine into single dataframe with all technical replicates
  combined <- bind_rows(reps_df)

  # Rename some protein groups
  if (length(relabel) > 0) {
    # Rename lazy eval function
    rename_lazy = interp(~str_replace_all(protein_desc, relabel), relabel = relabel)
    combined %>%
      mutate_(.dots = setNames(list(rename_lazy), c("protein_desc"))) -> combined
  }

  # Check if we should normalize to total area for a given replicate
  # This accounts for any variability in how the sample was injected
  message("Quantitating...")
  if (normalize == TRUE) {
    # Lazy mutate list
    normalize_lazy = list(~sum(area, na.rm = T), ~ area/total_area)
    combined %>%
      group_by_(~tech_rep) %>%
      mutate_(.dots = setNames(normalize_lazy, c("total_area", "area"))) %>%
      ungroup() -> combined
  }

  # Quantitate using top three most abundant areas
  combined %>%
    group_by_(~protein_desc) %>%
    do_(~merge_top_peptides(., num_reps, match_peps = match_peps)) -> combined

  return(combined)

}
