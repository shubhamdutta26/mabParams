test_that("calculate_mass throws error if not type character or illegal character", {
  expect_error(calculate_mass(1), "HC/ LC sequence input needs to be a character vector of length 1")
  expect_error(calculate_mass(TRUE), "HC/ LC sequence input needs to be a character vector of length 1")
  expect_error(calculate_mass(c("MRGM", "MGRM")), "HC/ LC sequence input needs to be a character vector of length 1")
  expect_error(calculate_mass("MRXM"), "Input contains:X; Only one letter amino acid codes should be in the input without spaces.")
  expect_error(calculate_mass("MRX M"), "Input contains:X,  ; Only one letter amino acid codes should be in the input without spaces.")
  expect_error(calculate_mass("MR}X M"), "Input contains:}, X,  ; Only one letter amino acid codes should be in the input without spaces.")
})

test_that("calculate_mass prints accurate reduced peptide masses daltons", {
  peptide <- "QVTLKESGPGILQPSQTLSLTCSFSGFSLRTSGMGV"
  outcome <- calculate_mass(peptide)
  expected <- 3715.22
  expect_equal(outcome, expected)
})

test_that("calculate_mass prints accurate reduced peptide masses daltons", {
  peptide <- "QVTLKESGPGILQPSQTLSLTCSFSGFSLRTSGMGVGWIRQPSGKGLEWLAHIWWDDDKRYNPALKSRLTISKDTSSNQVFLKIASVDTADTATYYCAQINPAWFAYWGQGTLVTVSSASTKGPSVFPLAPSSRSTSESTAALGCLVKDYFPEPVTVSWNSGSLTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYVCNVNHKPSNTKVDKRVEIKTCGGGSKPPTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSQEDPDVKFNWYVNGAEVHHAQTKPRETQYNSTYRVVSVLTVTHQDWLNGKEYTCKVSNKALPAPIQKTISKDKGQPREPQVYTLPPSREELTKNQVSLTCLVKGFYPSDIVVEWESSGQPENTYKTTPPVLDSDGSYFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSVSPGK"
  outcome <- calculate_mass(peptide)
  expected <- 49322.19
  expect_equal(outcome, expected)
})
