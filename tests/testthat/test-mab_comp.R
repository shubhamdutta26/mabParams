test_that("mab_comp throws error if not type character or illegal character", {
  expect_error(mab_comp(1, 1), "HC/ LC sequence input needs to be a character vector of length 1")
  expect_error(mab_comp(TRUE, FALSE), "HC/ LC sequence input needs to be a character vector of length 1")
  expect_error(mab_comp(c("MRGM", "MGRM"), "MGRM"), "HC/ LC sequence input needs to be a character vector of length 1")
  expect_error(mab_comp("MRXM", "MRGM"), "Input contains:X; Only one letter amino acid codes should be in the input without spaces.")
  expect_error(mab_comp("MRGM","MRX M"), "Input contains:X,  ; Only one letter amino acid codes should be in the input without spaces.")
  expect_error(mab_comp("MR}X M", "MR]X M"), "Input contains:}, X,  ; Only one letter amino acid codes should be in the input without spaces.")
  expect_error(mab_comp("MRGM"), 'argument "lc_seq" is missing, with no default')
  expect_error(mab_comp(lc_seq = "MRGM"), 'argument "hc_seq" is missing, with no default')
})

test_that("mab_comp returns the correct data", {
  hc = "QVTLKESGPGILQPSQTLSLTCSFSGFSLRTSGMGVGWIRQPSGKGLEWLAHIWWDDDKRYNPALKSRLTISKDTSSNQVFLKIASVDTADTATYYCAQINPAWFAYWGQGTLVTVSSASTKGPSVFPLAPSSRSTSESTAALGCLVKDYFPEPVTVSWNSGSLTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYVCNVNHKPSNTKVDKRVEIKTCGGGSKPPTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSQEDPDVKFNWYVNGAEVHHAQTKPRETQYNSTYRVVSVLTVTHQDWLNGKEYTCKVSNKALPAPIQKTISKDKGQPREPQVYTLPPSREELTKNQVSLTCLVKGFYPSDIVVEWESSGQPENTYKTTPPVLDSDGSYFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSVSPGK"
  lc = "DIVLTQSPASLAVSLGQRATISCKASQSVDFDGDSFMNWYQQKPGQPPKLLIYTTSNLESGIPARFSASGSGTDFTLNIHPVEEEDTATYYCQQSNEDPYTFGGGTKLELKRAVAAPSVFIFPPSEDQVKSGTVSVVCLLNNFYPREASVKWKVDGVLKTGNSQESVTEQDSKDNTYSLSSTLTLSSTDYQSHNVYACEVTHQGLSSPVTKSFNRGEC"
  expected <- readr::read_csv("test_data/mab_comp.csv", col_types = "ccddc", show_col_types = FALSE)
  expect_identical(mab_comp(hc, lc), expected)
})
