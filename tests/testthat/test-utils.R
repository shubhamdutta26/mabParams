test_that("check_input throws error if not type character or illegal character", {
  expect_error(check_input(1), "This needs to be a character vector of length 1")
  expect_error(check_input(TRUE), "This needs to be a character vector of length 1")
  expect_error(check_input(c("MRGM", "MGRM")), "This needs to be a character vector of length 1")
  expect_error(check_input("MRXM"), "Input contains:X; Only one letter amino acid codes should be in the input without spaces.")
  expect_error(check_input("MRX M"), "Input contains:X,  ; Only one letter amino acid codes should be in the input without spaces.")
  expect_error(check_input("MR}X M"), "Input contains:}, X,  ; Only one letter amino acid codes should be in the input without spaces.")
})

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

test_that("get_seq_element_comp returns tbl with element compositions", {
  input <- c("arg", "cys", "gly", "met")
  output <- tibble::tibble(molecule = c("arg", "cys", "gly", "met"),
                           carbon = c(6,3,2,5),
                           hydrogen = c(12,5,3,9),
                           nitrogen = c(4,1,1,1),
                           oxygen = c(1,1,1,1),
                           sulphur = c(0,1,0,1))
  result <- get_seq_element_comp(input)

  expect_identical(result, output)
})
