# arg parameters used for testing mab_params
hc_seq <- "QVTLKESGPGILQPSQTLSLTCSFSGFSLRTSGMGVGWIRQPSGKGLEWLAHIWWDDDKRYNPALKSRLTISKDTSSNQVFLKIASVDTADTATYYCAQINPAWFAYWGQGTLVTVSSASTKGPSVFPLAPSSRSTSESTAALGCLVKDYFPEPVTVSWNSGSLTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYVCNVNHKPSNTKVDKRVEIKTCGGGSKPPTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSQEDPDVKFNWYVNGAEVHHAQTKPRETQYNSTYRVVSVLTVTHQDWLNGKEYTCKVSNKALPAPIQKTISKDKGQPREPQVYTLPPSREELTKNQVSLTCLVKGFYPSDIVVEWESSGQPENTYKTTPPVLDSDGSYFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSVSPGK"
lc_seq <- "DIVLTQSPASLAVSLGQRATISCKASQSVDFDGDSFMNWYQQKPGQPPKLLIYTTSNLESGIPARFSASGSGTDFTLNIHPVEEEDTATYYCQQSNEDPYTFGGGTKLELKRAVAAPSVFIFPPSEDQVKSGTVSVVCLLNNFYPREASVKWKVDGVLKTGNSQESVTEQDSKDNTYSLSSTLTLSSTDYQSHNVYACEVTHQGLSSPVTKSFNRGEC"
mab <- "CD16"
chem_mod <- "C39H67N5O7"
cyclized <- TRUE
clipped <- TRUE

# Test begins
test_that("mab_params calculates correct masses for all FALSE and hc chemical modifications", {

  expected <- readr::read_csv("test_data/py_all_false_hc_mod.csv", show_col_types = FALSE)$mass
  actual <- round(mab_params(hc_seq, lc_seq, hc_chem_mod = chem_mod)$mass, digits = 0)

  expect_equal(actual, expected)
})

test_that("mab_params calculates correct masses for all FALSE and lc chemical modifications", {

  expected <- readr::read_csv("test_data/py_all_false_lc_mod.csv", show_col_types = FALSE)$mass
  actual <- round(mab_params(hc_seq, lc_seq, lc_chem_mod = chem_mod)$mass, digits = 0)

  expect_equal(actual, expected)
})

test_that("mab_params calculates correct masses for all FALSE and no chemical modifications", {

  expected <- readr::read_csv("test_data/py_all_false_no_mod.csv", show_col_types = FALSE)$mass
  actual <- round(mab_params(hc_seq, lc_seq)$mass, digits = 0)

  expect_equal(actual, expected)
})

test_that("mab_params calculates correct masses for all FALSE and both chemical modifications", {

  expected <- readr::read_csv("test_data/py_all_false_with_mod.csv", show_col_types = FALSE)$mass
  actual <- round(mab_params(hc_seq, lc_seq, hc_chem_mod = chem_mod, lc_chem_mod = chem_mod)$mass, digits = 0)

  expect_equal(actual, expected)
})

test_that("mab_params calculates correct masses for all TRUE and hc chemical modifications", {

  expected <- readr::read_csv("test_data/py_all_true_hc_mod.csv", show_col_types = FALSE)$mass
  actual <- round(mab_params(hc_seq, lc_seq, T, T, T, hc_chem_mod = chem_mod)$mass, digits = 0)

  expect_equal(actual, expected)
})

test_that("mab_params calculates correct masses for all TRUE and lc chemical modifications", {

  expected <- readr::read_csv("test_data/py_all_true_lc_mod.csv", show_col_types = FALSE)$mass
  actual <- round(mab_params(hc_seq, lc_seq, T, T, T, lc_chem_mod = chem_mod)$mass, digits = 0)

  expect_equal(actual, expected)
})

test_that("mab_params calculates correct masses for all TRUE and no chemical modifications", {

  expected <- readr::read_csv("test_data/py_all_true_no_mod.csv", show_col_types = FALSE)$mass
  actual <- round(mab_params(hc_seq, lc_seq, T, T, T)$mass, digits = 0)

  expect_equal(actual, expected)
})

test_that("mab_params calculates correct masses for all TRUE and both chemical modifications", {

  expected <- readr::read_csv("test_data/py_all_true_with_mod.csv", show_col_types = FALSE)$mass
  actual <- round(mab_params(hc_seq, lc_seq, T, T, T, hc_chem_mod = chem_mod, lc_chem_mod = chem_mod)$mass, digits = 0)

  expect_equal(actual, expected)
})

# Errors
test_that("mab_params throws error if not args are not valid input", {
  expect_error(mab_params(1,1),
               regexp = "chain sequence inputs must be a single")
  expect_error(mab_params(TRUE, FALSE),
               regexp = "chain sequence inputs must be a single")
  expect_error(mab_params(c("MRGM", "MGRM")),
               regexp = "chain sequence inputs must be a single")
  expect_error(mab_params("MRXM"),
               regexp = "Invalid character found in input sequence")
  expect_error(mab_params("MRX M"),
               regexp = "Sequence contains spaces")
  expect_error(mab_params("MR}XM"),
               regexp = "Invalid character found in input sequence")
  expect_error(mab_params("MRGM", "MRGM", 1,1,1),
               regexp = "must be logical")
  expect_error(mab_params("MRGM", "MRGM", hc_chem_mod = 1),
               "The heavy chain chemical modification needs to be a single character string.")
  expect_error(mab_params("MRGM", "MRGM", lc_chem_mod = 1),
               "The light chain chemical modification needs to be a single character string.")
  expect_error(mab_params("MRGM", "MRGM", hc_chem_mod = 1, lc_chem_mod = 1),
               "The heavy chain chemical modification needs to be a single character string.")
})
