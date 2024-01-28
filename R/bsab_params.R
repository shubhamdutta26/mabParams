basb_params <- function (hc_seq_1,
                         hc_seq_2,
                         lc_seq_1,
                         lc_seq_2,
                         hc_cyclized_1 = FALSE,
                         hc_cyclized_2 = FALSE,
                         hc_clipped_1 = FALSE,
                         hc_clipped_2 = FALSE,
                         lc_glycosylation = FALSE,
                         hc_chem_mod_1 = NA,
                         hc_chem_mod_2 = NA,
                         lc_chem_mod_1 = NA,
                         lc_chem_mod_2 = NA,
                         n_disulphides = 16L,
                         n_hc_disulphides = 4L,
                         n_lc_disulphides = 2L,
                         mab = "mab1") {
  # Read data-------------------------------------------------------------------
  element_composition <- readr::read_csv("inst/extdata/element_composition.csv",
                                         col_types = "ciiiii",
                                         col_names = TRUE)
  glycans <- readr::read_csv("inst/extdata/glycans.csv",
                             col_types = "ciiiiicc",
                             col_names = TRUE)
  element_symbol <- readr::read_csv("inst/extdata/element_symbols.csv",
                                    col_types = "cc",
                                    col_names = TRUE)

  # Check sequence and generate tibble with atomic composition------------------
  seq_list <- list(hc_seq_1, hc_seq_2, lc_seq_1, lc_seq_2)
  purrr::map(seq_list, check_input)
  boolean_args_list <- list(hc_cyclized_1, hc_cyclized_2, hc_clipped_1,
                            hc_clipped_2, lc_glycosylation)
  purrr::map(boolean_args_list, check_boolean_args)
}
