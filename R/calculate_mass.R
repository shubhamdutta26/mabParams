#' Calculate mass of peptides (reduced) in daltons
#'
#' @param seq One peptide sequence without spaces
#'
#' @return The numeric mass in daltons
#' @export
#'
#' @examples
#' sequence <- "QVTLKESGPGILQPSQTLSLTCSFSGFSLRTSGMGV"
#' mass <- calculate_mass(sequence)
#' print(mass)
calculate_mass <- function (seq) {
  check_input(seq)
  seq_three <- to_aa_three_letter(seq)
  seq_base <- c(seq_three, "water")
  seq_base_count_table <- count_molecules(seq_base)
  seq_to_elements <- get_seq_element_comp(seq_base)
  seq_elements_tibble <- seq_to_elements  %>%
    dplyr::left_join(seq_base_count_table, by = "molecule") %>%
    dplyr::mutate(dplyr::across(2:(ncol(.)-1), ~.x*n)) %>%
    dplyr::select(-n) %>%
    dplyr::select(2:ncol(.)) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(.x, na.rm = TRUE)))
  return(paste(round(calculate_mass_from_elements_tbl(seq_elements_tibble),
                     digits = 2), "Da", sep = " "))
}
