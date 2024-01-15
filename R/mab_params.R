#' Calculate mass parameters of monoclonal antibodies
#'
#' @param hc_seq A string with hc amino acids
#' @param lc_seq A string with lc amino acids
#' @param hc_cyclized A boolean
#' @param hc_clipped A boolean
#' @param lc_glycosylation A boolean
#' @param hc_chem_mod A string
#' @param lc_chem_mod A string
#' @param n_disulphides Number of total disulphides
#' @param n_hc_disulphides Number of total HC disulphides
#' @param n_lc_disulphides Number of total LC disulphides
#' @param mab Name of mab; default is mab1
#' @return vector
#' @export
#'
#' @examples
mab_params <- function (hc_seq,
                        lc_seq,
                        hc_cyclized = F,
                        hc_clipped = F,
                        lc_glycosylation = F,
                        hc_chem_mod = NA,
                        lc_chem_mod = NA,
                        n_disulphides = 16L,
                        n_hc_disulphides = 4L,
                        n_lc_disulphides = 2L,
                        mab = "mab1") {

  # HC MASS CALCULATION---------------------------------------------------------
  hc_mass_params <- calculate_hc_mass(hc_seq, hc_cyclized, hc_clipped,
                                      hc_chem_mod, n_hc_disulphides, mab)

  # LC MASS CALCULATION---------------------------------------------------------
  lc_mass_params <- calculate_lc_mass(lc_seq, lc_glycosylation, lc_chem_mod,
                                      n_lc_disulphides, mab)

  # INTACT MASS NO GLYCOFORMS---------------------------------------------------
  hc_three_1 <- to_aa_three_letter(hc_seq)
  hc_first_aa_1 <- hc_three_1[1]
  if (hc_cyclized == TRUE) {
    col_name_cyclized <- paste0("Pyro", ifelse(hc_first_aa_1 == "gln", "Q", "E"))
  } else {
    col_name_cyclized <- "cyclized"
  }

  cli_cyc_tbl <- hc_mass_params[[2]]
  lc_mod_elements_tibble <- lc_mass_params[[2]]
  n_hydrogen <- tibble::tibble(elements = "hydrogen",
                               number = n_disulphides * -2L)

  hc_count_tbl <- cli_cyc_tbl %>%
    dplyr::mutate(dplyr::across(2:(ncol(.)-1), ~.x*n)) %>%
    dplyr::select(-n) %>%
    dplyr::select(2:ncol(.)) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(.x, na.rm = TRUE))) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "elements", values_to = "hc") %>%
    dplyr::mutate(hc_twice = hc * 2) %>%
    dplyr::select(-hc)

  lc_count_tbl <- lc_mod_elements_tibble %>%
    dplyr::mutate(dplyr::across(2:(ncol(.)-1), ~.x*n)) %>%
    dplyr::select(-n) %>%
    dplyr::select(2:ncol(.)) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(.x, na.rm = TRUE))) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "elements", values_to = "lc") %>%
    dplyr::mutate(lc_twice = lc * 2) %>%
    dplyr::select(-lc)

  intact_mass <- purrr::reduce(list(hc_count_tbl, lc_count_tbl, n_hydrogen),
                                dplyr::full_join, by = "elements") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(total_count = sum(hc_twice, lc_twice, number, na.rm = TRUE)) %>%
    dplyr::select(elements, total_count) %>%
    tidyr::pivot_wider(names_from = elements, values_from = total_count) %>%
    dplyr::mutate(mass = calculate_mass_from_elements_tbl(.),
                  {{col_name_cyclized}} := ifelse(hc_cyclized == TRUE, "Yes", "No"),
                  clipping = ifelse(hc_clipped == TRUE, "Yes", "No")) %>%
    dplyr::mutate(antibody = mab, chain = NA, reduction = "intact",
                  modifications = paste(hc_chem_mod, lc_chem_mod, sep = "; "),
                  glycans = "No", .before = 1)

  # INTACT MASS WITH GLYCOFORMS-------------------------------------------------
  element_composition <- readr::read_csv("inst/extdata/element_composition.csv",
                                         col_types = "ciiiii",
                                         col_names = TRUE)

  glycans_1 <- readr::read_csv("inst/extdata/glycans.csv",
                             col_types = "ciiiiicc",
                             col_names = TRUE)
  hc <- glycans_1 %>%
    dplyr::filter(show_hc == "Y") %>%
    dplyr::select(-c(show_hc, show_lc))

  lc <- glycans_1 %>%
    dplyr::filter(show_lc == "Y") %>%
    dplyr::select(-c(show_hc, show_lc))

  combinations <- combine_and_sum(hc, lc, names(hc)[1], names(hc)[-1], "_hc", "_lc") %>%
    dplyr::mutate(hc_lc_modifications = paste0(glycan_name_hc, " + ", glycan_name_hc,
                                               " + ", glycan_name_lc, " + ",
                                               glycan_name_lc, sep = ""), .before = 1) %>%
    dplyr::select(hc_lc_modifications, dplyr::starts_with("sum_")) %>%
    dplyr::mutate(n = 2) %>%
    dplyr::mutate(dplyr::across(2:(ncol(.)-1), ~.x*n)) %>%
    select(-n) %>%
    dplyr::rename_with(~gsub("sum_", "", .x, fixed = TRUE)) %>%
    tidyr::pivot_longer(-hc_lc_modifications,
                        names_to = "molecule",
                        values_to = "n")

  split_combinations <- split(combinations, f = combinations$hc_lc_modifications)

  combinations_full <- split_combinations %>%
    purrr::map(dplyr::select, -1) %>%
    purrr::map(dplyr::left_join, element_composition, by = "molecule") %>%
    purrr::map(dplyr::relocate, n, .after = dplyr::last_col()) %>%
    # purrr::map(dplyr::bind_rows, cli_cyc_tbl, lc_mod_elements_tibble) %>%
    purrr::map(dplyr::mutate_if, is.numeric, ~.x*n) %>%
    purrr::map(dplyr::select, -dplyr::last_col()) %>%
    purrr::map(dplyr::select, 2:dplyr::last_col()) %>%
    purrr::map(dplyr::summarise_all, ~ sum(.x, na.rm = TRUE)) %>%
    purrr::map(tidyr::pivot_longer, dplyr::everything(), names_to = "elements", values_to = "glycans")

  intact_mass_glyco <- combinations_full %>%
    purrr::map(dplyr::full_join, hc_count_tbl, by = "elements") %>%
    purrr::map(dplyr::full_join, lc_count_tbl, by = "elements") %>%
    purrr::map(dplyr::full_join, n_hydrogen, by = "elements") %>%
    # purrr::map(dplyr::rowwise) %>%
    purrr::map(~dplyr::mutate(., total_sum = rowSums(.[,2:5], na.rm = TRUE))) %>%
    purrr::map(dplyr::select, c(elements, total_sum)) %>%
    purrr::map(tidyr::pivot_wider, names_from = elements, values_from = total_sum)

  final_glyco_intact <- dplyr::bind_rows(intact_mass_glyco) %>%
    dplyr::mutate(glycans = names(intact_mass_glyco), .before = 1) %>%
    dplyr::mutate(mass = calculate_mass_from_elements_tbl(.),
                  {{col_name_cyclized}} := ifelse(hc_cyclized == TRUE, "Yes", "No"),
                  clipping = ifelse(hc_clipped == TRUE, "Yes", "No")) %>%
    dplyr::mutate(antibody = mab, chain = NA, reduction = "intact",
                  modifications = paste(hc_chem_mod, lc_chem_mod, sep = "; "),
                  .before = 1)


  # RETURN TIBBLE
  return(dplyr::bind_rows(intact_mass, final_glyco_intact, hc_mass_params[[1]],lc_mass_params[[1]]))
}
