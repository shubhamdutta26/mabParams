test_that("parse_chem_formula returns error if first character is a number", {
  formula <- "1H23N2OS"
  expect_error(parse_chem_formula(formula), "The formula cannot begin with a number or contain special characters")
})

test_that("parse_chem_formula returns error if any spaces in formula", {
  formula <- "1H23 N2OS"
  expect_error(parse_chem_formula(formula), "The formula cannot begin with a number or contain special characters")
})

test_that("parse_chem_formula returns error if any special chars in formula", {
  formula <- "1H23(N2OS)2"
  expect_error(parse_chem_formula(formula), "The formula cannot begin with a number or contain special characters")
})

test_that("parse_chem_formula returns error if any special chars in formula - 2", {
  formula <- "1H23*N2OS2"
  expect_error(parse_chem_formula(formula), "The formula cannot begin with a number or contain special characters")
})

test_that("parse_chem_formula returns error if any special chars in formula - 3", {
  formula <- "1H23N2OS2/"
  expect_error(parse_chem_formula(formula), "The formula cannot begin with a number or contain special characters")
})

test_that("parse_chem_formula returns error if any special chars in formula - 4", {
  formula <- "$1H23N2OS2"
  expect_error(parse_chem_formula(formula), "The formula cannot begin with a number or contain special characters")
})

test_that("parse_chem_formula returns the correct element composition", {
  formula <- "CHNOS"
  parsed <- parse_chem_formula(formula)
  expected <- tibble::tibble(
    element = c("C", "H", "N", "O", "S"),
    count = 1
  )

  expect_identical(parsed, expected)
})

test_that("parse_chem_formula returns the correct element composition of C1H1N1O1S1", {
  formula <- "C1H1N1O1S1"
  parsed <- parse_chem_formula(formula)
  expected <- tibble::tibble(
    element = c("C", "H", "N", "O", "S"),
    count = 1
  )

  expect_identical(parsed, expected)
})

test_that("parse_chem_formula returns the correct element composition of C10H12N15O1S20", {
  formula <- "C10H12N15O1S20"
  parsed <- parse_chem_formula(formula)
  expected <- tibble::tibble(
    element = c("C", "H", "N", "O", "S"),
    count = c(10, 12, 15, 1, 20)
  )

  expect_identical(parsed, expected)
})

test_that("parse_chem_formula returns the correct element composition of C102H123N156O12S201", {
  formula <- "C102H123N156O12S201"
  parsed <- parse_chem_formula(formula)
  expected <- tibble::tibble(
    element = c("C", "H", "N", "O", "S"),
    count = c(102, 123, 156, 12, 201)
  )

  expect_identical(parsed, expected)
})

test_that("parse_chem_formula returns the correct element composition of C1025H1230N1561O123S2011", {
  formula <- "C1025H1230N1561O123S2011"
  parsed <- parse_chem_formula(formula)
  expected <- tibble::tibble(
    element = c("C", "H", "N", "O", "S"),
    count = c(1025, 1230, 1561, 123, 2011)
  )

  expect_identical(parsed, expected)
})

test_that("parse_chem_formula returns the correct element composition", {
  formula <- "CoHeNOSe"
  parsed <- parse_chem_formula(formula)
  expected <- tibble::tibble(
    element = c("Co", "He", "N", "O", "Se"),
    count = 1
  )

  expect_identical(parsed, expected)
})

test_that("parse_chem_formula returns the correct element composition of Co1He1N1O1Se1", {
  formula <- "Co1He1N1O1Se1"
  parsed <- parse_chem_formula(formula)
  expected <- tibble::tibble(
    element = c("Co", "He", "N", "O", "Se"),
    count = 1
  )

  expect_identical(parsed, expected)
})

test_that("parse_chem_formula returns the correct element composition of Co10He12N15O1Se20", {
  formula <- "Co10He12N15O1Se20"
  parsed <- parse_chem_formula(formula)
  expected <- tibble::tibble(
    element = c("Co", "He", "N", "O", "Se"),
    count = c(10, 12, 15, 1, 20)
  )

  expect_identical(parsed, expected)
})

test_that("parse_chem_formula returns the correct element composition of Co102He123N156O12Se201", {
  formula <- "Co102He123N156O12Se201"
  parsed <- parse_chem_formula(formula)
  expected <- tibble::tibble(
    element = c("Co", "He", "N", "O", "Se"),
    count = c(102, 123, 156, 12, 201)
  )

  expect_identical(parsed, expected)
})

test_that("parse_chem_formula returns the correct element composition of Co1025He1230N1561O123Se2011", {
  formula <- "Co1025He1230N1561O123Se2011"
  parsed <- parse_chem_formula(formula)
  expected <- tibble::tibble(
    element = c("Co", "He", "N", "O", "Se"),
    count = c(1025, 1230, 1561, 123, 2011)
  )

  expect_identical(parsed, expected)
})

