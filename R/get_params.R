#' XXX
#'
#' @param seq A string with amino acids
#'
#' @return vector
#' @export
#'
#' @examples
get_params <- function (seq) {
  check_input(seq)
  aa_three_vector <- to_aa_three_letter(seq)
  # return(aa_three_vector)
  count_table <- count_molecules(aa_three_vector)
  # return(count_table)
  element_table <- get_seq_element_comp(aa_three_vector)
  # return(element_table)
  joined <- element_table %>%
    dplyr::left_join(count_table, by = "molecule")
  # return(joined)
  get_total_mass(joined)
}
