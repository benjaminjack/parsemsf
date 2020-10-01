
#' Make a table of peptide areas
#'
#' Areas under each peptide peak that can be used downstream for quantitation.
#' See \code{\link{quantitate}} for protein quantitation.
#'
#' @inheritParams  make_pep_table
#'
#' @return A data frame containing peptide areas for peptides at or above the minimum confidence level.
#'
#' \item{peptide_id}{a unique peptide ID}
#' \item{spectrum_id}{a unique spectrum ID}
#' \item{protein_desc}{protein description from reference database used to assign peptides to protein groups, parsed according to \code{prot_regex}}
#' \item{sequence}{amino acid sequence (does not show post-translational modifications)}
#' \item{area}{area under peptide peak}
#' \item{mass}{peptide mass}
#' \item{m_z}{mass-to-charge ratio}
#' \item{charge}{peptide charge}
#' \item{intensity}{peak intensity; useful if no area is available}
#' \item{first_scan}{first scan in which peptide appears}
#'
#' @export
#'
#' @examples
#' make_area_table(parsemsf_example("test_db.msf"))
make_area_table <- function(msf_file, min_conf = "High", prot_regex = "^>([a-zA-Z0-9._]+)\\b", collapse = TRUE) {

  # Access MSF database file
  # my_db <- src_sqlite(msf_file)
  my_db <- DBI::dbConnect(RSQLite::SQLite(), msf_file)
  on.exit(DBI::dbDisconnect(my_db))

  Events <- tbl(my_db, "Events") %>%
    select(-"RT", -"LeftRT", -"RightRT", -"SN", -"FileID", -"Intensity")

  EventAreaAnnotations <- tbl(my_db, "EventAreaAnnotations") %>%
    select("EventID", "QuanResultID")

  PrecursorIonAreaSearchSpectra <- tbl(my_db, "PrecursorIonAreaSearchSpectra")

  # Here is where all the scan information, including mass and charge
  SpectrumHeaders <- tbl(my_db, "SpectrumHeaders") %>%
    select("SpectrumID", "FirstScan", "Mass", "Charge", "RetentionTime", "UniqueSpectrumID") %>%
    collect(n = Inf)

  # Grab intensities
  MassPeaks <- tbl(my_db, "MassPeaks") %>%
    select("MassPeakID", "Intensity") %>%
    collect(n = Inf)

  # Here are all the areas
  # m/z is really just one of many m/z's that could be retrieved from the msf file.
  # All of the m/z's are calculated slightly differently. I just picked one since
  # I don't use m/z's for my work.
  events_joined <- Events %>%
    inner_join(EventAreaAnnotations, by = "EventID") %>%
    inner_join(PrecursorIonAreaSearchSpectra, by = "QuanResultID") %>%
    collect(n = Inf)  %>%
    group_by(.data$SearchSpectrumID) %>%
    summarize(
      Area = sum(.data$Area),
      m_z = mean(.data$Mass)
    )

  # NOTE: mass and m/z are probably not correct right now! I have check them in more detail

  # Join areas to the spectrum IDs. Not all spectrum IDs have areas!
  # I haven't figured out what the difference between UniqueSpectrumID and SpectrumID are,
  # except that MassPeakID's correspond to UniqueSpectrumID ONLY.
  spectra <- left_join(SpectrumHeaders, events_joined, by = c("SpectrumID" = "SearchSpectrumID")) %>%
    inner_join(MassPeaks, by = c("UniqueSpectrumID" = "MassPeakID"))

  pep_table <- make_pep_table(msf_file, min_conf, prot_regex, collapse)

  auc_table <- spectra %>%
    # Join peptide info to mass/area/charge/etc.
    right_join(pep_table, by = c("SpectrumID" = "spectrum_id")) %>%
    # Rename columns for more consistent column naming
    select(
      "peptide_id" = "peptide_id",
      "spectrum_id" = "SpectrumID",
      "protein_desc" = "protein_desc",
      "sequence" = "sequence",
      "area" = "Area",
      "mass" = "Mass",
      "m_z" = "m_z",
      "charge" = "Charge",
      "intensity" = "Intensity",
      "first_scan" = "FirstScan"
    )

  return(auc_table)
}
