
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
    select(SpectrumID, FirstScan, Mass, Charge, RetentionTime) %>%
    collect()
  
  # Grab intensities
  MassPeaks <- tbl(my_db, "MassPeaks") %>%
    select(MassPeakID, Intensity) %>%
    collect()
  
  # Here are all the areas
  events_joined <- inner_join(Events, EventAreaAnnotations) %>% 
    inner_join(PrecursorIonAreaSearchSpectra) %>%
    collect()  %>%
    rename(m_z = Mass)
  
  # NOTE: mass and m/z are probably not correct right now! I have check them in more detail
  
  # Join areas to the spectrum IDs. Not all spectrum IDs have areas!
  spectra <- left_join(SpectrumHeaders, events_joined, by = c("SpectrumID" = "SearchSpectrumID")) %>%
    inner_join(MassPeaks, by = c("SpectrumID" = "MassPeakID"))
  
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
           FirstScan)
  
  return(auc_table)
  
}