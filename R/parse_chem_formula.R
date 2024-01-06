element_symbol <- c(
  bromine = "Br",
  calcium = "Ca",
  carbon = "C",
  chlorine = "Cl",
  fluorine = "F",
  hydrogen = "H",
  iodine = "I",
  lithium = "Li",
  nitrogen = "N",
  oxygen = "O",
  phosphorus = "P",
  potassium = "K",
  selenium = "Se",
  sodium = "Na",
  sulphur = "S"
)

parse_chem_formula <- function(formula,
                               element = "element",
                               count = "count",
                               total_count = "total_count") {
  #formula <- "NH3"
  if (stringr::str_detect(formula, pattern = " ") == TRUE ||
      stringr::str_detect(formula, pattern = "[$&+,:;=?@#|'<>.^*()%!-]") == TRUE) {
    stop("Cannot have spaces or special characters, Only A-Z, a-z, and 0-9 allowed")
  }
  # Split the formula into individual characters
  chars <- stringr::str_split(formula, "", simplify = TRUE)

  # Initialize vectors to store the elements and their counts
  symbols <- c()
  counts <- c()

  # Initialize a variable to store the current element
  current_element <- ""

  # Loop over each character in the formula
  for (i in 1:length(chars)) {
    # If the character is an uppercase letter, it's the start of a new element
    if (grepl("[A-Z]", chars[i])) {
      # If there's a current element, add it to the symbols vector with a count of 1
      if (current_element != "") {
        symbols <- c(symbols, current_element)
        counts <- c(counts, 1)
      }

      # Set the current element to the new element
      current_element <- chars[i]
    } else if (grepl("[a-z]", chars[i])) {
      # If the character is a lowercase letter, it's part of the current element
      current_element <- paste0(current_element, chars[i])
    } else if (grepl("[0-9]", chars[i])) {
      # If the character is a number, it's the count of the current element
      symbols <- c(symbols, current_element)
      counts <- c(counts, as.integer(chars[i]))

      # Reset the current element
      current_element <- ""
    }
  }

  # If there's a current element at the end of the formula, add it to the symbols vector with a count of 1
  if (current_element != "") {
    symbols <- c(symbols, current_element)
    counts <- c(counts, 1)
  }

  # Create a dataframe from the vectors
  df <- tibble::tibble(!!element := symbols, !!count := counts) %>%
    dplyr::group_by(element) %>%
    dplyr::summarise(!!total_count := sum(count)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from = element, values_from = total_count) %>%
    dplyr::mutate(molecule = formula, .before = 1) %>%
    dplyr::rename(dplyr::any_of(element_symbol))

  # Return the dataframe
  return(df)
}
