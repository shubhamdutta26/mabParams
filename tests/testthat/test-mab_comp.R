test_that("mab_comp throws error if not type character or illegal character", {
  expect_error(mab_comp(1, 1),
               regexp = "chain sequence inputs must be a single")
  expect_error(mab_comp(TRUE, FALSE),
               regexp = "chain sequence inputs must be a single")
  expect_error(mab_comp(c("MRGM", "MGRM"),
                        regexp = "chain sequence inputs must be a single"))
  expect_error(mab_comp("MRXM", "MRGM"),
               regexp = "Invalid character found in input sequence")
  expect_error(mab_comp("MRGM","MRX M"),
               regexp = "Sequence contains spaces")
  expect_error(mab_comp("MR}XM", "MR]XM"),
               regexp = "Invalid character found in input sequence")
  expect_error(mab_comp("MRGM"), 'argument "lc_seq" is missing, with no default')
  expect_error(mab_comp(lc_seq = "MRGM"), 'argument "hc_seq" is missing, with no default')
})

test_that("mab_comp returns the correct data", {
  hc = "QVTLKESGPGILQPSQTLSLTCSFSGFSLRTSGMGVGWIRQPSGKGLEWLAHIWWDDDKRYNPALKSRLTISKDTSSNQVFLKIASVDTADTATYYCAQINPAWFAYWGQGTLVTVSSASTKGPSVFPLAPSSRSTSESTAALGCLVKDYFPEPVTVSWNSGSLTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYVCNVNHKPSNTKVDKRVEIKTCGGGSKPPTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSQEDPDVKFNWYVNGAEVHHAQTKPRETQYNSTYRVVSVLTVTHQDWLNGKEYTCKVSNKALPAPIQKTISKDKGQPREPQVYTLPPSREELTKNQVSLTCLVKGFYPSDIVVEWESSGQPENTYKTTPPVLDSDGSYFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSVSPGK"
  lc = "DIVLTQSPASLAVSLGQRATISCKASQSVDFDGDSFMNWYQQKPGQPPKLLIYTTSNLESGIPARFSASGSGTDFTLNIHPVEEEDTATYYCQQSNEDPYTFGGGTKLELKRAVAAPSVFIFPPSEDQVKSGTVSVVCLLNNFYPREASVKWKVDGVLKTGNSQESVTEQDSKDNTYSLSSTLTLSSTDYQSHNVYACEVTHQGLSSPVTKSFNRGEC"
  expected <- readr::read_csv("test_data/mab_comp.csv", col_types = "ccddc", show_col_types = FALSE)
  expect_identical(mab_comp(hc, lc), expected)
})
