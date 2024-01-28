#' Converts chemical formula into a tibble with names and element counts. It does
#' not check if the letters match any chemical formula. For e.g. if X is present
#' in the formula. which is not in the periodic table it doesnot throw an error.
#' There cannot be any spaces or special characters in the formula.
#'
#' @param formula A character vector og length 1 (e.g. C2H6O for ethanol; not C2H5OH)
#'
#' @return A tibble/ dataframe
#' @export
#'
#' @examples
#' formula <- "C6H12O6" # glucose
#' parse_chem_formula(formula)
parse_chem_formula <- function(formula) {
  # Check if first character is a number or any special characters
  if (grepl("^[0-9]$", substr(formula, 1, 1)) == TRUE || grepl("[^a-zA-Z0-9]", formula) == TRUE ) {
    stop("The formula cannot begin with a number or contain special characters.", call. = FALSE)
  }
  # Use regular expression to match element symbols and their counts
  matches <- regmatches(formula, gregexpr("[A-Z][a-z]?\\d*()", formula))[[1]]

  # Initialize vectors to store elements and counts
  elements <- gsub("\\d", "", matches)
  counts <- gsub("\\D", "", matches)
  counts[counts == ""] <- "1"
  counts <- as.numeric(counts)

  # Construct and return dataframe
  return(tibble::tibble(element = elements, count = counts))
}
