library(parsemsf)
context("Make peptide table")

pep_col_names <- c("PeptideID",
                   "SpectrumID",
                   "Sequence",
                   "Proteins",
                   "ProteinID",
                   "PEP",
                   "q-Value")

test_that("make_pep_table creates a data frame with the correct column names", {
  expect_equal(colnames(make_pep_table(parsemsf_example("test_db.msf"))), pep_col_names)
})
