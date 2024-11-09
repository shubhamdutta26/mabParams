# check_input() tests
test_that("check_input throws error if not type character or illegal character", {
  expect_error(check_input(1),
               regexp = "chain sequence inputs must be a single")
  expect_error(check_input(TRUE),
               regexp = "chain sequence inputs must be a single")
  expect_error(check_input(c("MRGM", "MGRM")),
               regexp = "chain sequence inputs must be a single")
  expect_error(check_input("MRXM"),
               regexp = "Invalid character found in input sequence")
  expect_error(check_input("MRX M"),
               regexp = "Sequence contains spaces")
  expect_error(check_input("MR}M"),
               regexp = "Invalid character found in input sequence")
  expect_silent(check_input("EVQLVESGGGLVQPGGSLRLSCAASGFTFS"))
  expect_silent(check_input("MRRM"))
  expect_silent(check_input("ACDEFF"))
})

# check_boolean_args() tests
test_that("check_boolean_args accepts valid boolean inputs", {
  # Should not throw any errors
  expect_no_error(check_boolean_args(TRUE, TRUE, TRUE))
  expect_no_error(check_boolean_args(FALSE, FALSE, FALSE))
  expect_no_error(check_boolean_args(TRUE, FALSE, TRUE))
})

test_that("check_boolean_args accepts NULL values", {
  # Test all combinations with NULL
  expect_no_error(check_boolean_args(NULL, NULL, NULL))
  expect_no_error(check_boolean_args(TRUE, NULL, NULL))
  expect_no_error(check_boolean_args(NULL, TRUE, NULL))
  expect_no_error(check_boolean_args(NULL, NULL, TRUE))
  expect_no_error(check_boolean_args(TRUE, TRUE, NULL))
  expect_no_error(check_boolean_args(TRUE, NULL, TRUE))
  expect_no_error(check_boolean_args(NULL, TRUE, TRUE))
})

test_that("check_boolean_args rejects non-boolean values", {
  # Test numeric values
  expect_error(
    check_boolean_args(1, TRUE, TRUE),
    "The following argument\\(s\\) must be logical \\(TRUE/FALSE\\): cyclized"
  )

  # Test character values
  expect_error(
    check_boolean_args(TRUE, "FALSE", TRUE),
    "The following argument\\(s\\) must be logical \\(TRUE/FALSE\\): clipped"
  )

  # Test other R objects
  expect_error(
    check_boolean_args(TRUE, TRUE, list()),
    "The following argument\\(s\\) must be logical \\(TRUE/FALSE\\): glycosylation"
  )
})

test_that("check_boolean_args correctly identifies multiple invalid arguments", {
  # Test multiple invalid arguments
  expect_error(
    check_boolean_args(1, "FALSE", TRUE),
    "The following argument\\(s\\) must be logical \\(TRUE/FALSE\\): cyclized, clipped"
  )

  expect_error(
    check_boolean_args("TRUE", 2, "FALSE"),
    "The following argument\\(s\\) must be logical \\(TRUE/FALSE\\): cyclized, clipped, glycosylation"
  )
})

# to_aa_three_letter() tests
test_that("to_aa_three_letter returns a character vector of input string upper or lower case)", {
  expect_identical(to_aa_three_letter("MRGMGGC"), c("met", "arg", "gly", "met", "gly", "gly", "cys"))
  expect_identical(to_aa_three_letter("MrGMggC"), c("met", "arg", "gly", "met", "gly", "gly", "cys"))
})

test_that("count_molecules returns a class list and counts each element in the vector", {
  input <- c("met", "arg", "gly", "met", "gly", "gly", "cys")
  output <- tibble::tibble(molecule = c("arg", "cys", "gly", "met"),
                           n = c(1L, 1L, 3L, 2L))
  expect_identical(count_molecules(input), output)
  expect_identical(typeof(count_molecules(input)), "list")
})

# get_seq_element_comp() tests
test_that("get_seq_element_comp returns tbl with element compositions", {
  input <- c("arg", "cys", "gly", "met")
  output <- tibble::tibble(molecule = c("arg", "cys", "gly", "met"),
                           carbon = c(6,3,2,5),
                           hydrogen = c(12,5,3,9),
                           nitrogen = c(4,1,1,1),
                           oxygen = c(1,1,1,1),
                           sulphur = c(0,1,0,1))
  result <- get_seq_element_comp(input)
  expect_equal(result, output)
})

# Create a minimal test dataset
test_element_composition <- data.frame(
  molecule = c("water", "carbon", "ala", "gly"),
  carbon = c(0L, 1L, 3L, 2L),
  hydrogen = c(2L, 0L, 5L, 3L),
  nitrogen = c(0L, 0L, 1L, 1L),
  oxygen = c(1L, 0L, 1L, 1L),
  sulphur = c(0L, 0L, 0L, 0L)
)

test_that("get_seq_element_comp handles basic functionality", {
  # Test single molecule lookup
  result <- get_seq_element_comp("water")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(result$molecule, "water")
  expect_equal(result$carbon, 0)
  expect_equal(result$hydrogen, 2)

  # Test multiple molecule lookup
  result <- get_seq_element_comp(c("ala", "gly"))
  expect_equal(nrow(result), 2)
  expect_equal(result$molecule, c("ala", "gly"))
  expect_equal(result$carbon, c(3, 2))
})

test_that("get_seq_element_comp preserves data frame attributes", {
  result <- get_seq_element_comp("water")

  # Check column names
  expect_equal(
    colnames(result),
    c("molecule", "carbon", "hydrogen", "nitrogen", "oxygen", "sulphur")
  )

  # Check column types
  expect_type(result$molecule, "character")
  expect_type(result$carbon, "integer")
  expect_type(result$hydrogen, "integer")
  expect_type(result$nitrogen, "integer")
  expect_type(result$oxygen, "integer")
  expect_type(result$sulphur, "integer")
})

# calculate_mass_from_elements_tbl() tests
test_that("calculate_mass_from_elements_tbl calculates correct masses", {
  # Test case 1: Simple molecule (Water - H2O)
  water <- data.frame(
    hydrogen = c(2),
    oxygen = c(1)
  )
  expected_mass_water <- 2 * 1.00794075 + 15.99940492
  expect_equal(
    calculate_mass_from_elements_tbl(water),
    expected_mass_water,
    tolerance = 1e-7
  )

  # Test case 2: Multiple compounds (Methanol and Ethanol)
  alcohols <- data.frame(
    carbon = c(1, 2),
    hydrogen = c(4, 6),
    oxygen = c(1, 1)
  )
  expected_masses <- c(
    (1 * 12.0107359 + 4 * 1.00794075 + 1 * 15.99940492),  # Methanol
    (2 * 12.0107359 + 6 * 1.00794075 + 1 * 15.99940492)   # Ethanol
  )
  expect_equal(
    calculate_mass_from_elements_tbl(alcohols),
    expected_masses,
    tolerance = 1e-7
  )

  # Test case 3: Handles NA values
  compounds_with_na <- data.frame(
    carbon = c(1, NA),
    hydrogen = c(4, 6),
    oxygen = c(1, 1)
  )
  expect_false(any(is.na(calculate_mass_from_elements_tbl(compounds_with_na))))

  # Test case 5: Empty data frame
  empty_df <- data.frame(
    carbon = numeric(0),
    hydrogen = numeric(0)
  )
  expect_equal(
    length(calculate_mass_from_elements_tbl(empty_df)),
    0
  )
})
