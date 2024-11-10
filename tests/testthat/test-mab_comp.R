library(dplyr, quietly = TRUE, warn.conflicts = FALSE)

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
})

test_that("mab_comp handles counts correctly", {
  # Test sequence with known counts
  hc <- "QQQ"  # 3 glutamines
  lc <- "QQ"   # 2 glutamines

  result <- mab_comp(hc, lc)

  # Check counts for each chain type
  heavy_gln <- result %>%
    dplyr::filter(chain == "heavy", amino_acid == "gln") %>%
    dplyr::pull(n)
  expect_equal(heavy_gln, 3)  # Single heavy chain

  light_gln <- result %>%
    dplyr::filter(chain == "light", amino_acid == "gln") %>%
    dplyr::pull(n)
  expect_equal(light_gln, 2)  # Single light chain

  both_gln <- result %>%
    dplyr::filter(chain == "both", amino_acid == "gln") %>%
    dplyr::pull(n)
  expect_equal(both_gln, 10)  # (3+2) * 2 for complete antibody
})

test_that("mab_comp handles percentages correctly", {
  hc <- "QVQLVQ"  # 6 residues
  lc <- "DIQMTQ"  # 6 residues

  result <- mab_comp(hc, lc) %>%
    group_by(chain) %>%
    summarise(total_perc = sum(aa_percent)) %>%
    ungroup()

  expect_equal(result$chain, c("both", "heavy", "light"))
  expect_equal(result$total_perc, c(1,1,1), tolerance = 0.01)
})

test_that("mab_comp handles invalid input appropriately", {
  expect_error(mab_comp("ABC123", "DEF456"))  # Invalid amino acids
  expect_error(mab_comp("", "QVQLVQ"))        # Empty string
  expect_error(mab_comp(NULL, "QVQLVQ"))      # NULL input
  expect_error(mab_comp("QVQLVQ"))
})

test_that("mab_comp properly joins with aa_type data", {
  hc <- "QVQLVQ"
  lc <- "DIQMTQ"

  result <- mab_comp(hc, lc)

  # Check that all amino acids have polarity and charge information
  expect_true(all(!is.na(result$polarity)))
  expect_true(all(!is.na(result$net_charge_7.4)))

  # Check specific classifications
  asp_row <- result[result$amino_acid == "asp", ]
  expect_equal(asp_row$net_charge_7.4, c("Negative", "Negative"))
  expect_equal(asp_row$polarity, c("Brønsted base", "Brønsted base"))
})
