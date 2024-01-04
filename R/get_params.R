#' XXX
#'
#' @param hc_seq A string with hc amino acids
#' @param hc_cyclized Boolean
#' @param hc_clipped Boolean
#' @return vector
#' @export
#'
#' @examples
get_params <- function (hc_seq, hc_cyclized, hc_clipped) {
  check_input(hc_seq)
  aa_three_vector <- to_aa_three_letter(hc_seq)
  # return(aa_three_vector)
  count_table <- count_molecules(aa_three_vector)
  # return(count_table)
  element_table <- get_seq_element_comp(aa_three_vector)
  # return(element_table)
  joined <- element_table %>%
    dplyr::left_join(count_table, by = "molecule")
  # return(joined)
  n = "n"
  counts = "counts"
  mass = "mass"
  each_mass = "each_mass"
  hc_mass <- calculate_mass(joined, n, counts, mass, each_mass)

  # Cyclization check
  symbol = "symbol"
  mass = "mass"
  nh3_mass <- element_mass %>%
    dplyr::filter(symbol == "N") %>%
    dplyr::pull(mass) +
    (element_mass %>%
    dplyr::filter(symbol == "H") %>%
    dplyr::pull(mass))*3

  h2o_mass <- element_mass %>%
    dplyr::filter(symbol == "H") %>%
    dplyr::pull(mass) +
    (element_mass %>%
       dplyr::filter(symbol == "O") %>%
       dplyr::pull(mass))*2

  if (hc_cyclized == TRUE) {
    if (stringr::str_starts(hc_seq, "Q") == TRUE) {
      hc_mass = hc_mass - nh3_mass
    } else if (stringr::str_starts(hc_seq, "E") == TRUE) {
      hc_mass = hc_mass - h2o_mass
    } else {
      stop(paste0("HC can only be cyclized if it starts with Q (PyroQ) or E (PyroE)"),
           call. = FALSE)
    }
  }

  # Lysine clipping check
  K_mass = (6*12.0107359)+(12*1.00794075)+(2*14.00670321)+(1*15.99940492)

  if (hc_clipped == TRUE) {
    if (stringr::str_ends(hc_seq, "K") == TRUE) {
      hc_mass = hc_mass - K_mass
    }
  }

  return(hc_mass)





}
