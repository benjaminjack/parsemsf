library(parsemsf)
context("Map peptides")

pep_col_names <- c("PeptideID",
                   "SpectrumID",
                   "Pep_seq",
                   "Proteins",
                   "ProteinID",
                   "PEP",
                   "q-Value",
                   "Prot_seq",
                   "start",
                   "end")

test_that("map_peptides creates a data frame with the correct column names", {
  expect_equal(colnames(map_peptides(parsemsf_example("test_db.msf"))), pep_col_names)
})
