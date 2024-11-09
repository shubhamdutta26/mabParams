#' Standard one-letter amino acid codes
#' @keywords internal
# aa_one_letter <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
#                    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

#' @keywords internal
#'
#' Internal function to validate antibody HC/LC sequence input
#'
#' Validates whether the input sequence for antibody heavy chain (HC) or light chain (LC)
#' is in the correct format and contains only valid one-letter amino acid codes.
#'
#' @param one_letter_input A character string containing the amino acid sequence
#' @return NULL invisibly. Raises an error if validation fails
check_input <- function(one_letter_input) {
  # Input type and length validation
  if (!is.character(one_letter_input) || length(one_letter_input) != 1) {
    rlang::abort(
      "Heavy and light chain sequence inputs must be a single character string each.",
      call = NULL
    )
  }

  # Check for spaces first
  if (grepl("\\s", one_letter_input)) {
    rlang::abort(
      "Sequence contains spaces. Please remove all spaces from the input.",
      call = NULL
    )
  }

  # Convert to uppercase and split into individual characters
  # Fix: Using strsplit instead of str_split to ensure proper character vector
  split_aa <- strsplit(toupper(one_letter_input), "")[[1]]

  # Check for invalid amino acids using the standard aa_one_letter vector
  invalid_aa <- split_aa[!split_aa %in% aa_one_letter]

  if (length(invalid_aa) > 0) {
    invalid_list <- paste(invalid_aa, collapse = ", ")
    rlang::abort(
      c(
        "Invalid character found in input sequence.",
        "x" = sprintf(
          "Invalid amino acid: %s",
          invalid_list
        ),
        "i" = "Only standard one-letter amino acid codes (A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V) are allowed."
      ),
      call = NULL
    )
  }

  invisible(NULL)
}

#' @keywords internal
#'
#' Check Boolean Arguments
#'
#' @description
#' Validates that specified arguments (Cyclization, clipping, and glycosylation)
#' are boolean values (TRUE/FALSE). Provides specific error messages for any
#' non-boolean arguments.
#'
#' @param cyclization A logical value indicating whether cyclization is enabled
#' @param clipping A logical value indicating whether clipping is enabled
#' @param glycosylation A logical value indicating whether glycosylation is enabled
#'
#' @return Nothing. Called for side effects.
#' @throws error if any argument is not boolean
#' @noexport
check_boolean_args <- function(cyclized = NULL, clipped = NULL, glycosylation = NULL) {
  # Create named list of arguments for checking
  args_to_check <- list(
    cyclized = cyclized,
    clipped = clipped,
    glycosylation = glycosylation
  )

  # Check each non-NULL argument and collect non-boolean ones
  invalid_args <- names(args_to_check)[sapply(args_to_check, function(x) {
    !is.null(x) && !is.logical(x)
  })]

  # If any invalid arguments found, construct and throw error message
  if (length(invalid_args) > 0) {
    error_msg <- sprintf(
      "The following argument(s) must be logical (TRUE/FALSE): %s",
      paste(invalid_args, collapse = ", ")
    )
    rlang::abort(error_msg, call = NULL)
  }
}

#' Convert One-Letter to Three-Letter Amino Acid Codes
#'
#' @description
#' Converts amino acid sequences from one-letter codes to their corresponding
#' three-letter codes. Handles both uppercase and lowercase inputs.
#'
#' @param one_letter_input A character string containing one-letter amino acid codes
#'   (e.g., "ACGT" or "acgt")
#'
#' @return A character vector of three-letter amino acid codes
#'
#' @details
#' The function accepts a string of one-letter amino acid codes and converts each
#' letter to its corresponding three-letter code. It is case-insensitive and will
#' handle both upper and lowercase inputs. Invalid amino acid codes will result in
#' an error.
#'
#' @examples
#' \dontrun{
#' to_aa_three_letter("ACD")  # Returns c("ala", "cys", "asp")
#' to_aa_three_letter("acd")  # Same result
#' }
#'
#' @keywords internal
to_aa_three_letter <- function(one_letter_input) {

  # Split string into individual characters and convert to uppercase
  aa_chars <- strsplit(one_letter_input, "")[[1]] |>
    toupper()
  # Convert to three-letter codes
  unname(aa_named[aa_chars])
}

#' Count occurrences of each unique molecule
#'
#' @description
#' Creates a tibble containing counts of how many times each unique molecule appears
#' in the input vector.
#'
#' @param molecule A vector containing molecule identifiers
#'
#' @return A tibble with two columns:
#'   \itemize{
#'     \item molecule: The unique molecule values
#'     \item n: Count of how many times each molecule appears
#'   }
#'
#' @examples
#' molecules <- c("ATP", "GTP", "ATP", "CTP")
#' count_molecules(molecules)
#'
#' @keywords internal
count_molecules <- function(molecule) {
  tibble::as_tibble(table(molecule))
}

#' Get elemental composition for specified molecules
#'
#' Internal function that retrieves the elemental composition data for a given
#' set of molecules from the element_composition dataset.
#'
#' @param seq Character vector of molecule names to look up
#'   \itemize{
#'     \item molecule: Character. Name of the molecule
#'     \item carbon: Integer. Number of carbon atoms
#'     \item hydrogen: Integer. Number of hydrogen atoms
#'     \item nitrogen: Integer. Number of nitrogen atoms
#'     \item oxygen: Integer. Number of oxygen atoms
#'     \item sulphur: Integer. Number of sulfur atoms
#'   }
#'
#' @return A data frame containing rows from element_composition matching the
#'   input molecules
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Get composition for multiple molecules
#' get_seq_element_comp(c("ala", "gly", "ser"))
#'
#' # Get composition for a single molecule
#' get_seq_element_comp("water")
#' }
get_seq_element_comp <- function(seq) {
  element_composition[element_composition$molecule %in% seq, ]
}

check_chem_mod_args <- function(hc_chem_mod = NA, lc_chem_mod = NA) {
  if (!is.na(hc_chem_mod)) {
    if (typeof(hc_chem_mod) != "character" || length(hc_chem_mod) > 1) {
      stop("The chemical modifications needs to be a character vector of length 1.", call. = FALSE)
    }
    if (grepl("^[0-9]$", substr(hc_chem_mod, 1, 1)) == TRUE || grepl("[^a-zA-Z0-9]", hc_chem_mod) == TRUE) {
      stop("The chemical modifications cannot begin with a number or contain special characters.", call. = FALSE)
    }
  }
  if (!is.na(lc_chem_mod)) {
    if (typeof(lc_chem_mod) != "character" || length(lc_chem_mod) > 1) {
      stop("The chemical modifications needs to be a character vector of length 1.", call. = FALSE)
    }
    if (grepl("^[0-9]$", substr(lc_chem_mod, 1, 1)) == TRUE || grepl("[^a-zA-Z0-9]", lc_chem_mod) == TRUE) {
      stop("The chemical modifications cannot begin with a number or contain special characters.", call. = FALSE)
    }
  }
}

#' Calculate Molecular Mass from Element Composition
#'
#' @description
#' Takes a data frame containing element quantities and calculates the total molecular mass
#' for each row by multiplying element quantities with their respective atomic masses
#' and summing the results.
#'
#' @param df A data frame where column names match element names (e.g., "carbon", "hydrogen")
#'           and values represent the number of atoms of each element
#'
#' @return A numeric vector containing the calculated molecular mass for each row
#'
#' @details
#' The function uses a predefined `element_mass_list` containing atomic masses for elements.
#' Only columns matching element names in this list are used in calculations.
#' Non-matching columns are ignored. Missing values (NA) are treated as 0 in the sum.
#'
#' @examples
#' \dontrun{
#' compounds <- data.frame(
#'   carbon = c(1, 2),
#'   hydrogen = c(4, 6),
#'   oxygen = c(1, 1)
#' )
#' calculate_mass_from_elements_tbl(compounds)
#' }
#'
#' @keywords internal
calculate_mass_from_elements_tbl <- function(df) {
  # Multiply each column by the corresponding element mass
  element_cols <- intersect(names(df), names(element_mass_list))
  # Apply the element mass multiplication
  df[element_cols] <- lapply(element_cols, function(col) {
    df[[col]] * element_mass_list[[col]]
  })
  # Calculate row-wise mass sum across numeric columns
  df$mass_full <- rowSums(df[sapply(df, is.numeric)], na.rm = TRUE)
  # Return the calculated mass
  df$mass_full
}

combine_and_sum <- function(df1, df2, char_col, num_cols, suffix_1, suffix_2) {
  comb <- expand.grid(df1[[char_col]], df2[[char_col]])
  names(comb) <- c(paste0(char_col, suffix_1), paste0(char_col, suffix_2))

  for (num_col in num_cols) {
    comb[[paste0(num_col, suffix_1)]] <- df1[[num_col]][match(comb[[paste0(char_col, suffix_1)]], df1[[char_col]])]
    comb[[paste0(num_col, suffix_2)]] <- df2[[num_col]][match(comb[[paste0(char_col, suffix_2)]], df2[[char_col]])]
    comb[[paste0("sum_", num_col)]] <- comb[[paste0(num_col, suffix_1)]] + comb[[paste0(num_col, suffix_2)]]
  }

  return(comb)
}
