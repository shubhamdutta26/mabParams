calculate_hc_mass <- function (seq,
                               cyclized,
                               clipped,
                               chem_mod = NA,
                               n_hc_disulphides = 4L,
                               mab) {
  # seq = "QVTLKESGPGILQPSQTLSLTCSFSGFSLRTSGMGVGWIRQPSGKGLEWLAHIWWDDDKRYNPALKSRLTISKDTSSNQVFLKIASVDTADTATYYCAQINPAWFAYWGQGTLVTVSSASTKGPSVFPLAPSSRSTSESTAALGCLVKDYFPEPVTVSWNSGSLTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYVCNVNHKPSNTKVDKRVEIKTCGGGSKPPTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSQEDPDVKFNWYVNGAEVHHAQTKPRETQYNSTYRVVSVLTVTHQDWLNGKEYTCKVSNKALPAPIQKTISKDKGQPREPQVYTLPPSREELTKNQVSLTCLVKGFYPSDIVVEWESSGQPENTYKTTPPVLDSDGSYFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSVSPGK"
  # mab = "mab1"
  # chem_mod = "C39H67N5O7"

  # Read data-------------------------------------------------------------------
  element_composition <- readr::read_csv("inst/extdata/element_composition.csv",
                                         col_types = "ciiiii",
                                         col_names = TRUE)
  glycans <- readr::read_csv("inst/extdata/glycans.csv",
                             col_types = "ciiiiicc",
                             col_names = TRUE)
  element_symbol <- c(bromine = "Br", calcium = "Ca", carbon = "C",
                      chlorine = "Cl", fluorine = "F", hydrogen = "H",
                      iodine = "I", lithium = "Li", nitrogen = "N",
                      oxygen = "O", phosphorus = "P", potassium = "K",
                      selenium = "Se", sodium = "Na", sulphur = "S")

  # Check sequence and generate tibble with atomic composition------------------
  check_input(seq)
  hc_three <- to_aa_three_letter(seq)
  hc_base <- c(hc_three, "water")
  hc_base_count_table <- count_molecules(hc_base)
  hc_to_elements <- get_seq_element_comp(hc_base)
  hc_elements_tibble <- hc_to_elements  %>%
    dplyr::left_join(hc_base_count_table, by = "molecule")

  # Check mods------------------------------------------------------------------
  if (!is.na(chem_mod)) {
    parsed_chem_tbl <- parse_chem_formula(chem_mod) %>%
      tidyr::pivot_wider(names_from = element, values_from = count) %>%
      dplyr::mutate(molecule = chem_mod, n = 1L) %>%
      dplyr::rename(dplyr::any_of(element_symbol))
  } else {
    parsed_chem_tbl <- NULL
  }

  hc_mod_elements_tibble <- hc_elements_tibble %>%
    dplyr::bind_rows(parsed_chem_tbl) %>%
    dplyr::mutate_all(~replace(., is.na(.), 0L))

  # Cyclization check-----------------------------------------------------------
  hc_first_aa <- hc_three[1]
  if (cyclized == TRUE & !(hc_first_aa %in% c("gln", "glu"))) {
    stop(paste0("HC can only be cyclized if it starts with Q (PyroQ) or E (PyroE)"),
         call. = FALSE)
  }

  # Clipping check--------------------------------------------------------------
  hc_last_aa <- hc_three[length(hc_three)]
  if (clipped == TRUE & hc_last_aa != "lys") {
    stop(paste0("HC can only be clipped if ends with Lysine (K)"),
         call. = FALSE)
  }

  # Cyclization and clipping calculations of HC full reduced MW-----------------
  hc_first_aa <- hc_three[1]
  hc_last_aa <- hc_three[length(hc_three)]
  if (cyclized == TRUE) {
    col_name_cyclized <- paste0("Pyro", ifelse(hc_first_aa == "gln", "Q", "E"))
  } else {
    col_name_cyclized <- "cyclized"
  }
  # cyclized = TRUE
  # clipped = TRUE

  if (cyclized == TRUE & clipped == TRUE) {
    if (hc_first_aa == "gln") {
      cli_cyc_tbl <- hc_mod_elements_tibble %>%
        dplyr::bind_rows(
          tibble::tibble(
            element_composition %>% dplyr::filter(molecule %in% c("lys", "nh3")),
            n = -1L
          )
        )
    } else {
      cli_cyc_tbl <- hc_mod_elements_tibble %>%
        dplyr::bind_rows(
          tibble::tibble(
            element_composition %>% dplyr::filter(molecule %in% c("lys", "water")),
            n = -1L
          )
        )
    }
  } else if (cyclized == TRUE & clipped == FALSE) {
    if (hc_first_aa == "gln") {
      cli_cyc_tbl <- hc_mod_elements_tibble %>%
        dplyr::bind_rows(
          tibble::tibble(
            element_composition %>% dplyr::filter(molecule == "nh3"),
            n = -1L
          )
        )
    } else {
      cli_cyc_tbl <- hc_mod_elements_tibble %>%
        dplyr::bind_rows(
          tibble::tibble(
            element_composition %>% dplyr::filter(molecule == "water"),
            n = -1L
          )
        )
    }
  } else if (cyclized == FALSE & clipped == TRUE) {
    cli_cyc_tbl <- hc_mod_elements_tibble %>%
      dplyr::bind_rows(
        tibble::tibble(
          element_composition %>% dplyr::filter(molecule == "lys"),
          n = -1L
        )
      )
  } else if (cyclized == FALSE & clipped == FALSE) {
    cli_cyc_tbl <- hc_mod_elements_tibble
  }

  hc_mab_full_reduced <- cli_cyc_tbl %>%
    dplyr::mutate(dplyr::across(2:(ncol(.)-1), ~.x*n)) %>%
    dplyr::select(-n) %>%
    dplyr::select(2:ncol(.)) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(.x, na.rm = TRUE))) %>%
    dplyr::mutate(mass = calculate_mass_from_elements_tbl(.),
                  {{col_name_cyclized}} := ifelse(cyclized == TRUE, "Yes", "No"),
                  clipping = ifelse(clipped == TRUE, "Yes", "No")) %>%
    dplyr::mutate(antibody = mab, chain = "heavy", reduction = "full",
                  modifications = chem_mod, glycans = "No", .before = 1)

  # Calculations for partially reduced HC---------------------------------------
  if (cyclized == TRUE & clipped == TRUE) {
    if (hc_first_aa == "gln") {
      partial_red_tbl <- hc_mod_elements_tibble %>%
        dplyr::bind_rows(
          tibble::tibble(
            element_composition %>% dplyr::filter(molecule %in% c("lys", "nh3")),
            n = -1L
          )
        ) %>%
        dplyr::bind_rows(
          tibble::tibble(
            element_composition %>% dplyr::filter(molecule == "hydrogen"),
            n = -(2L * n_hc_disulphides)
          )
        )
    } else {
      partial_red_tbl <- hc_mod_elements_tibble %>%
        dplyr::bind_rows(
          tibble::tibble(
            element_composition %>% dplyr::filter(molecule %in% c("lys", "water")),
            n = -1L
          )
        ) %>%
        dplyr::bind_rows(
          tibble::tibble(
            element_composition %>% dplyr::filter(molecule == "hydrogen"),
            n = -(2L * n_hc_disulphides)
          )
        )
    }
  } else if (cyclized == TRUE & clipped == FALSE) {
    if (hc_first_aa == "gln") {
      partial_red_tbl <- hc_mod_elements_tibble %>%
        dplyr::bind_rows(
          tibble::tibble(
            element_composition %>% dplyr::filter(molecule == "nh3"),
            n = -1L
          )
        ) %>%
        dplyr::bind_rows(
          tibble::tibble(
            element_composition %>% dplyr::filter(molecule == "hydrogen"),
            n = -(2L * n_hc_disulphides)
          )
        )
    } else {
      partial_red_tbl <- hc_mod_elements_tibble %>%
        dplyr::bind_rows(
          tibble::tibble(
            element_composition %>% dplyr::filter(molecule == "water"),
            n = -1L
          )
        ) %>%
        dplyr::bind_rows(
          tibble::tibble(
            element_composition %>% dplyr::filter(molecule == "hydrogen"),
            n = -(2L * n_hc_disulphides)
          )
        )
    }
  } else if (cyclized == FALSE & clipped == TRUE) {
    partial_red_tbl <- hc_mod_elements_tibble %>%
      dplyr::bind_rows(
        tibble::tibble(
          element_composition %>% dplyr::filter(molecule == "lys"),
          n = -1L
        )
      ) %>%
      dplyr::bind_rows(
        tibble::tibble(
          element_composition %>% dplyr::filter(molecule == "hydrogen"),
          n = -(2L * n_hc_disulphides)
        )
      )
  } else if (cyclized == FALSE & clipped == FALSE) {
    partial_red_tbl <- hc_mod_elements_tibble %>%
      dplyr::bind_rows(
        tibble::tibble(
          element_composition %>% dplyr::filter(molecule == "hydrogen"),
          n = -(2L * n_hc_disulphides)
        )
      )
  }

  hc_mab_partial_reduced <- partial_red_tbl %>%
    dplyr::mutate(dplyr::across(2:(ncol(.)-1), ~.x*n)) %>%
    dplyr::select(-n) %>%
    dplyr::select(2:ncol(.)) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(.x, na.rm = TRUE))) %>%
    dplyr::mutate(mass = calculate_mass_from_elements_tbl(.),
                  {{col_name_cyclized}} := ifelse(cyclized == TRUE, "Yes", "No"),
                  clipping = ifelse(clipped == TRUE, "Yes", "No")) %>%
    dplyr::mutate(antibody = mab, chain = "heavy", reduction = "partial",
                  modifications = chem_mod, glycans = "No", .before = 1)

  # Add glycoforms to total reduced hc mass-------------------------------------
  hc_glycan_tbl <- glycans %>%
    dplyr::filter(show_hc == "Y") %>%
    dplyr::select(-c(show_hc, show_lc)) %>%
    tidyr::pivot_longer(-glycan_name,
                        names_to = "molecule",
                        values_to = "n")

  hc_glycan_individual <- split(hc_glycan_tbl, f = hc_glycan_tbl$glycan_name)

  hc_glycan_full <- hc_glycan_individual %>%
    purrr::map(dplyr::select, -1) %>%
    purrr::map(dplyr::left_join, element_composition, by = "molecule") %>%
    purrr::map(dplyr::relocate, n, .after = dplyr::last_col()) %>%
    purrr::map(dplyr::bind_rows, cli_cyc_tbl) %>%
    purrr::map(dplyr::mutate_if, is.numeric, ~.x*n) %>%
    purrr::map(dplyr::select, -dplyr::last_col()) %>%
    purrr::map(dplyr::select, 2:dplyr::last_col()) %>%
    purrr::map(dplyr::summarise_all, ~ sum(.x, na.rm = TRUE))

  hc_glycan_full_calculation <- dplyr::bind_rows(hc_glycan_full) %>%
    dplyr::mutate(glycans = names(hc_glycan_full), .before = 1) %>%
    dplyr::mutate(mass = calculate_mass_from_elements_tbl(.),
                  {{col_name_cyclized}} := ifelse(cyclized == TRUE, "Yes", "No"),
                  clipping = ifelse(clipped == TRUE, "Yes", "No"),
                  antibody = mab, chain = "heavy", reduction = "full",
                  modifications = chem_mod)

  hc_glycan_partial <- hc_glycan_individual %>%
    purrr::map(dplyr::select, -1) %>%
    purrr::map(dplyr::left_join, element_composition, by = "molecule") %>%
    purrr::map(dplyr::relocate, n, .after = dplyr::last_col()) %>%
    purrr::map(dplyr::bind_rows, partial_red_tbl) %>%
    purrr::map(dplyr::mutate_if, is.numeric, ~.x*n) %>%
    purrr::map(dplyr::select, -dplyr::last_col()) %>%
    purrr::map(dplyr::select, 2:dplyr::last_col()) %>%
    purrr::map(dplyr::summarise_all, ~ sum(.x, na.rm = TRUE))

  hc_glycan_partial_calculation <- dplyr::bind_rows(hc_glycan_partial) %>%
    dplyr::mutate(glycans = names(hc_glycan_partial), .before = 1) %>%
    dplyr::mutate(mass = calculate_mass_from_elements_tbl(.),
                  {{col_name_cyclized}} := ifelse(cyclized == TRUE, "Yes", "No"),
                  clipping = ifelse(clipped == TRUE, "Yes", "No"),
                  antibody = mab, chain = "heavy", reduction = "partial",
                  modifications = chem_mod)

 # Return parameters as tibble-------------------------------------------------
  hc_mass_params <- dplyr::bind_rows(hc_mab_full_reduced,
                                     hc_glycan_full_calculation,
                                     hc_mab_partial_reduced,
                                     hc_glycan_partial_calculation)

  return(list(hc_mass_params, cli_cyc_tbl))
}
