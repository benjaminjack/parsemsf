
make_auc_table <- function(msf_file, min_conf = "High", as_df = TRUE) {
  
  # Access MSF database file
  my_db <- src_sqlite(msf_file)
  
  Events <- tbl(my_db, "Events") %>%
    select(-RT, -LeftRT, -RightRT, -SN, -FileID, -Intensity)
  
  EventAreaAnnotations <- tbl(my_db, "EventAreaAnnotations") %>%
    select(EventID, QuanResultID)
  
  PrecursorIonAreaSearchSpectra <- tbl(my_db, "PrecursorIonAreaSearchSpectra")
  
  # Here is where all the scan information, including mass and charge
  SpectrumHeaders <- tbl(my_db, "SpectrumHeaders") %>%
    select(SpectrumID, FirstScan, Mass, Charge, RetentionTime, UniqueSpectrumID) %>%
    collect()
  
  # Grab intensities
  MassPeaks <- tbl(my_db, "MassPeaks") %>%
    select(MassPeakID, Intensity) %>%
    collect()
  
  # Here are all the areas
  # m/z is really just one of many m/z's that could be retrieved from the msf file.
  # All of the m/z's are calculated slightly differently. I just picked one since
  # I don't use m/z's for my work.
  events_joined <- inner_join(Events, EventAreaAnnotations) %>% 
    inner_join(PrecursorIonAreaSearchSpectra) %>%
    collect()  %>%
    group_by(SearchSpectrumID) %>%
    summarize(Area = sum(Area), m_z = mean(Mass))
  
  # NOTE: mass and m/z are probably not correct right now! I have check them in more detail
  
  # Join areas to the spectrum IDs. Not all spectrum IDs have areas!
  # I haven't figured out what the difference between UniqueSpectrumID and SpectrumID are,
  # except that MassPeakID's correspond to UniqueSpectrumID ONLY.
  spectra <- left_join(SpectrumHeaders, events_joined, by = c("SpectrumID" = "SearchSpectrumID")) %>%
    inner_join(MassPeaks, by = c("UniqueSpectrumID" = "MassPeakID"))
  
  pep_table <- make_pep_table(msf_file, min_conf)
  
  # Join peptide info to mass/area/charge/etc.
  auc_table <- right_join(spectra, pep_table, by=c("SpectrumID" = "SpectrumID")) %>% 
    select(PeptideID, 
           Sequence, 
           Proteins, 
           Area, 
           Mass, 
           m_z, 
           Charge, 
           Intensity, 
           FirstScan,
           SpectrumID)
  
  return(auc_table)
  
}