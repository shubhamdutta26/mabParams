#' Converts chemical formula into a tibble with names and element counts. Only
#' works for crtain elements. Add more elements if needed in the named vector
#' inside the function
#'
#' @param formula A character vector og length 1 (e.g. C2H5OH for ethanol)
#'
#' @return A tibble/ dataframe
#' @export
#'
#' @examples
parse_chem_formula <- function(formula) {
  # Use regular expression to match element symbols and their counts
  matches <- strsplit(formula, "(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)", perl=TRUE)[[1]]

  # Initialize vectors to store elements and counts
  elements <- matches[seq(1, length(matches), by=2)]
  counts <- as.numeric(matches[seq(2, length(matches), by=2)])

  # Construct dataframe
  df <- tibble::tibble(element = elements, count = counts)

  return(df)
}
