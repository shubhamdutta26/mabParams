#' Calculate mass parameters of monoclonal antibodies
#'
#' @param hc_seq A string containing heavy chain amino acids (one letter code)
#' @param lc_seq A string containing light chain amino acids (one letter code)
#' @param hc_cyclized A boolean value for heavy chain cyclization
#' @param hc_clipped A boolean value for heavy chain clipping
#' @param lc_glycosylation A boolean value for light chain glycosylation
#' @param hc_chem_mod A character vector of chemical formula of length 1
#' @param lc_chem_mod A character vector of chemical formula of length 1
#' @param n_disulphides Number of total disulphides
#' @param n_hc_disulphides Number of total heavy chain disulphides
#' @param n_lc_disulphides Number of total light chain disulphides
#' @param glycans A csv file that contains glycan parameters
#' @param mab Name of the antibody; default is "mab1"
#' @return A dataframe/ tibble with antibody mass parameters
#' @export
#'
#' @examples
#' heavy <- "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYA"
#' light <- "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSR"
#' glycans_path <- system.file("extdata", "glycans.csv", package = "mabParams")
#' params_tbl <- mab_params(heavy, light, glycans = glycans_path)
#' head(params_tbl)
mab_params <- function (hc_seq,
                        lc_seq,
                        hc_cyclized = FALSE,
                        hc_clipped = FALSE,
                        lc_glycosylation = FALSE,
                        hc_chem_mod = NA,
                        lc_chem_mod = NA,
                        n_disulphides = 16L,
                        n_hc_disulphides = 4L,
                        n_lc_disulphides = 2L,
                        glycans = NA,
                        mab = "mab1") {

  # Check chemical mod arg
  check_chem_mod_args(hc_chem_mod, lc_chem_mod)

  # Glycan input
  if (is.na(glycans)) {
    glycans_path <- system.file("extdata", "glycans.csv", package = "mabParams")
    glycans <- readr::read_csv(glycans_path,
                               col_types = "ciiiiicc",
                               col_names = TRUE)
  } else {
    glycans <- readr::read_csv(glycans,
                               col_types = "ciiiiicc",
                               col_names = TRUE)
  }

  # HC MASS CALCULATION---------------------------------------------------------
  hc_mass_params <- calculate_hc_mass(hc_seq, hc_cyclized, hc_clipped,
                                      hc_chem_mod, n_hc_disulphides, glycans, mab)

  # LC MASS CALCULATION---------------------------------------------------------
  lc_mass_params <- calculate_lc_mass(lc_seq, lc_glycosylation, lc_chem_mod,
                                      n_lc_disulphides, glycans, mab)

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
  # element_composition <- readr::read_csv("inst/extdata/element_composition.csv",
  #                                        col_types = "ciiiii",
  #                                        col_names = TRUE)
  hc <- glycans %>%
    dplyr::filter(show_hc == "Y") %>%
    dplyr::select(-c(show_hc, show_lc))


  lc_1 <- glycans %>%
    dplyr::filter(show_lc == "Y") %>%
    dplyr::select(-c(show_hc, show_lc))

  if (lc_glycosylation) {
    lc <- lc_1
  } else {
    lc <- lc_1 %>%
      dplyr::mutate(glycan_name = NA, HexNAc = 0, Hex = 0, dHex = 0, NeuAc = 0, NeuGc = 0)
  }

  combinations <- combine_and_sum(hc, lc, names(hc)[1], names(hc)[-1], "_hc", "_lc") %>%
    dplyr::mutate(hc_lc_modifications = paste0(glycan_name_hc, " + ", glycan_name_hc,
                                               " + ", glycan_name_lc, " + ",
                                               glycan_name_lc, sep = ""), .before = 1) %>%
    dplyr::select(hc_lc_modifications, dplyr::starts_with("sum_")) %>%
    dplyr::mutate(n = ifelse(lc_glycosylation == FALSE, 1, 2)) %>%
    dplyr::mutate(dplyr::across(2:(ncol(.)-1), ~.x*n)) %>%
    dplyr::select(-n) %>%
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
