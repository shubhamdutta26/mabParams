test_that("calculate_mass throws error if not type character or illegal character", {
  expect_error(calculate_mass(1),
               regexp = "chain sequence inputs must be a single")
  expect_error(calculate_mass(TRUE),
               regexp = "chain sequence inputs must be a single")
  expect_error(calculate_mass(c("MRGM", "MGRM")),
               regexp = "chain sequence inputs must be a single")
  expect_error(calculate_mass("MRXM"),
               regexp = "Invalid character found in input sequence")
  expect_error(calculate_mass("MRX M"),
               regexp = "Sequence contains spaces")
  expect_error(calculate_mass("MR}XM"),
               regexp = "Invalid character found in input sequence")
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
