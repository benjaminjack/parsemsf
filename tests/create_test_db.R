# "Proteins"
# "PeptidesProteins"
# "Events"
# "EventAreaAnnotations"
# "PrecursorIonAreaSearchSpectra"
# "SpectrumHeaders"
# "MassPeaks"
# "ProteinAnnotations"

library(dplyr)

prot_id = "807657"

my_db <- src_sqlite("../test_data/11-42_1min_5b.msf")

test_db <- src_sqlite("inst/extdata/test_db2.msf", create = T)

# Make SchemaInfo
SchemaInfo <- tbl(my_db, "SchemaInfo") %>% collect()
copy_to(test_db, SchemaInfo, temporary = FALSE)

# Make Proteins
Proteins <- tbl(my_db, "Proteins") %>% filter(ProteinID == prot_id) %>% collect()
copy_to(test_db, Proteins, temporary = FALSE)

# Make PeptidesProteins
PeptidesProteins <- tbl(my_db, "PeptidesProteins") %>% filter(ProteinID == prot_id) %>% collect()
copy_to(test_db, PeptidesProteins, temporary = FALSE)

# Make Peptides
pep_id_list <- PeptidesProteins$PeptideID
Peptides <- tbl(my_db, "Peptides") %>% filter(PeptideID %in% pep_id_list) %>% collect()
copy_to(test_db, Peptides, temporary = FALSE)

# Make ProteinAnnotations
ProteinAnnotations <- tbl(my_db, "ProteinAnnotations") %>% filter(ProteinID == prot_id) %>% collect()
copy_to(test_db, ProteinAnnotations, temporary = FALSE)

# Make CustomPeptides
CustomDataPeptides <- tbl(my_db, sql("SELECT FieldID, PeptideID, CAST(FieldValue as REAL) AS FieldValue FROM CustomDataPeptides")) %>%
  filter(PeptideID %in% pep_id_list) %>% collect()
copy_to(test_db, CustomDataPeptides, temporary = FALSE)

# Make CustomDataFields
CustomDataFields <- tbl(my_db, "CustomDataFields") %>% filter(FieldID %in% CustomDataPeptides$FieldID) %>% collect()
copy_to(test_db, CustomDataFields, temporary = FALSE)

# Make PrecursorIonAreaSearchSpectra
PrecursorIonAreaSearchSpectra <- tbl(my_db, "PrecursorIonAreaSearchSpectra") %>% filter(SearchSpectrumID %in% Peptides$SpectrumID) %>% collect()
copy_to(test_db, PrecursorIonAreaSearchSpectra, temporary = FALSE)

# Make EventAreaAnnotations
EventAreaAnnotations <- tbl(my_db, "EventAreaAnnotations") %>% filter(QuanResultID %in% PrecursorIonAreaSearchSpectra$QuanResultID) %>% collect()
copy_to(test_db, EventAreaAnnotations, temporary = FALSE)

# Make Events table
Events <- tbl(my_db, "Events") %>%
  filter(EventID %in% EventAreaAnnotations$EventID) %>%
  collect()
copy_to(test_db, Events, temporary = FALSE)

# Make SpectrumHeaders
SpectrumHeaders <- tbl(my_db, "SpectrumHeaders") %>% filter(SpectrumID %in% Peptides$SpectrumID) %>% collect()
copy_to(test_db, SpectrumHeaders, temporary = FALSE)

# Make MassPeaks
MassPeaks <- tbl(my_db, "MassPeaks") %>% filter(MassPeakID %in% SpectrumHeaders$UniqueSpectrumID) %>% collect()
copy_to(test_db, MassPeaks, temporary = FALSE)







