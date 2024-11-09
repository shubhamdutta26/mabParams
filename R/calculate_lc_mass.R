calculate_lc_mass <- function (seq,
                               glycosylation = FALSE,
                               chem_mod = NA,
                               n_lc_disulphides = 2L,
                               glycans,
                               mab) {

  # Read data-------------------------------------------------------------------
  # element_composition <- readr::read_csv("inst/extdata/element_composition.csv",
  #                                        col_types = "ciiiii",
  #                                        col_names = TRUE)
  # glycans <- readr::read_csv("inst/extdata/glycans.csv",
  #                            col_types = "ciiiiicc",
  #                            col_names = TRUE)
  # element_symbol <- readr::read_csv("inst/extdata/element_symbols.csv",
  #                                   col_types = "cc",
  #                                   col_names = TRUE)

  # seq = "DIVLTQSPASLAVSLGQRATISCKASQSVDFDGDSFMNWYQQKPGQPPKLLIYTTSNLESGIPARFSASGSGTDFTLNIHPVEEEDTATYYCQQSNEDPYTFGGGTKLELKRAVAAPSVFIFPPSEDQVKSGTVSVVCLLNNFYPREASVKWKVDGVLKTGNSQESVTEQDSKDNTYSLSSTLTLSSTDYQSHNVYACEVTHQGLSSPVTKSFNRGEC"
  # Check sequence and generate tibble with atomic composition------------------
  check_input(seq)
  check_boolean_args(glycosylation = glycosylation)
  lc_three <- to_aa_three_letter(seq)
  lc_base <- c(lc_three, "water")
  lc_base_count_table <- count_molecules(lc_base)
  lc_to_elements <- get_seq_element_comp(lc_base)
  lc_elements_tibble <- lc_to_elements  %>%
    dplyr::left_join(lc_base_count_table, by = "molecule")

  # Check mods------------------------------------------------------------------
  if (!is.na(chem_mod)) {
    parsed_chem_tbl <- parse_chem_formula(chem_mod) %>%
      dplyr::left_join(element_symbol, by = c("element" = "symbol")) %>%
      dplyr::select(-element) %>%
      dplyr::rename(element = element.y) %>%
      tidyr::pivot_wider(names_from = element, values_from = count) %>%
      dplyr::mutate(molecule = chem_mod, n = 1L)
  } else {
    parsed_chem_tbl <- NULL
  }

  lc_mod_elements_tibble <- lc_elements_tibble %>%
    dplyr::bind_rows(parsed_chem_tbl) %>%
    dplyr::mutate_all(~replace(., is.na(.), 0L))

  # Full reduction LC-----------------------------------------------------------
  lc_mab_full_reduced <- lc_mod_elements_tibble %>%
    dplyr::mutate(dplyr::across(2:(ncol(.)-1), ~.x*n)) %>%
    dplyr::select(-n) %>%
    dplyr::select(2:ncol(.)) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(.x, na.rm = TRUE))) %>%
    dplyr::mutate(mass = calculate_mass_from_elements_tbl(.)) %>%
    dplyr::mutate(antibody = mab, chain = "light", reduction = "full",
                  modifications = chem_mod, glycans = "No", .before = 1)

  # Partial reduction LC--------------------------------------------------------
  partial_red_tbl_lc <- lc_mod_elements_tibble %>%
    dplyr::bind_rows(
      tibble::tibble(
        element_composition %>% dplyr::filter(molecule == "hydrogen"),
        n = -(2L * n_lc_disulphides)
      )
    )

  lc_mab_partial_reduced <- partial_red_tbl_lc %>%
    dplyr::mutate(dplyr::across(2:(ncol(.)-1), ~.x*n)) %>%
    dplyr::select(-n) %>%
    dplyr::select(2:ncol(.)) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(.x, na.rm = TRUE))) %>%
    dplyr::mutate(mass = calculate_mass_from_elements_tbl(.)) %>%
    dplyr::mutate(antibody = mab, chain = "light", reduction = "partial",
                  modifications = chem_mod, glycans = "No", .before = 1)

  # Add glycoforms to total reduced lc mass-------------------------------------
  if (glycosylation == T) {
    lc_glycan_tbl <- glycans %>%
      dplyr::filter(show_lc == "Y") %>%
      dplyr::select(-c(show_hc, show_lc)) %>%
      tidyr::pivot_longer(-glycan_name,
                          names_to = "molecule",
                          values_to = "n")

    lc_glycan_individual <- split(lc_glycan_tbl, f = lc_glycan_tbl$glycan_name)

    lc_glycan_full <- lc_glycan_individual %>%
      purrr::map(dplyr::select, -1) %>%
      purrr::map(dplyr::left_join, element_composition, by = "molecule") %>%
      purrr::map(dplyr::relocate, n, .after = dplyr::last_col()) %>%
      purrr::map(dplyr::bind_rows, lc_mod_elements_tibble) %>%
      purrr::map(dplyr::mutate_if, is.numeric, ~.x*n) %>%
      purrr::map(dplyr::select, -dplyr::last_col()) %>%
      purrr::map(dplyr::select, 2:dplyr::last_col()) %>%
      purrr::map(dplyr::summarise_all, ~ sum(.x, na.rm = TRUE))

    lc_glycan_full_calculation <- dplyr::bind_rows(lc_glycan_full) %>%
      dplyr::mutate(glycans = names(lc_glycan_full), .before = 1) %>%
      dplyr::mutate(mass = calculate_mass_from_elements_tbl(.),
                    antibody = mab, chain = "light", reduction = "full",
                    modifications = chem_mod)

    lc_glycan_partial <- lc_glycan_individual %>%
      purrr::map(dplyr::select, -1) %>%
      purrr::map(dplyr::left_join, element_composition, by = "molecule") %>%
      purrr::map(dplyr::relocate, n, .after = dplyr::last_col()) %>%
      purrr::map(dplyr::bind_rows, partial_red_tbl_lc) %>%
      purrr::map(dplyr::mutate_if, is.numeric, ~.x*n) %>%
      purrr::map(dplyr::select, -dplyr::last_col()) %>%
      purrr::map(dplyr::select, 2:dplyr::last_col()) %>%
      purrr::map(dplyr::summarise_all, ~ sum(.x, na.rm = TRUE))

    lc_glycan_partial_calculation <- dplyr::bind_rows(lc_glycan_partial) %>%
      dplyr::mutate(glycans = names(lc_glycan_partial), .before = 1) %>%
      dplyr::mutate(mass = calculate_mass_from_elements_tbl(.),
                    antibody = mab, chain = "light", reduction = "partial",
                    modifications = chem_mod)
  } else {
    lc_glycan_full_calculation <- NULL
    lc_glycan_partial_calculation <- NULL
  }

  # Return parameters as tibble
  lc_mass_params <- dplyr::bind_rows(lc_mab_full_reduced,
                                     lc_glycan_full_calculation,
                                     lc_mab_partial_reduced,
                                     lc_glycan_partial_calculation)

  return(list(lc_mass_params, lc_mod_elements_tibble))

}
