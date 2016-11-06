library(parsemsf)
context("Make area table")

pep_col_names <- c("PeptideID",
                   "Sequence",
                   "Proteins",
                   "Area",
                   "Mass",
                   "m_z",
                   "Charge",
                   "Intensity",
                   "FirstScan",
                   "SpectrumID")

test_that("make_area_table creates a data frame with the correct column names", {
  expect_equal(colnames(make_area_table(parsemsf_example("test_db.msf"))), pep_col_names)
})
