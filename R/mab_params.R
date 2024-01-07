#' XXX
#'
#' @param hc_seq A string with hc amino acids
#' @param lc_seq A string with lc amino acids
#' @param hc_cyclized A boolean
#' @param hc_clipped A boolean
#' @param lc_glycosylation A boolean
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
                        hc_cyclized,
                        hc_clipped,
                        lc_glycosylation,
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
  # HC Mass calculation---------------------------------------------------------
  #hc_seq = "QVTLKESGPGILQPSQTLSLTCSFSGFSLRTSGMGVGWIRQPSGKGLEWLAHIWWDDDKRYNPALKSRLTISKDTSSNQVFLKIASVDTADTATYYCAQINPAWFAYWGQGTLVTVSSASTKGPSVFPLAPSSRSTSESTAALGCLVKDYFPEPVTVSWNSGSLTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYVCNVNHKPSNTKVDKRVEIKTCGGGSKPPTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSQEDPDVKFNWYVNGAEVHHAQTKPRETQYNSTYRVVSVLTVTHQDWLNGKEYTCKVSNKALPAPIQKTISKDKGQPREPQVYTLPPSREELTKNQVSLTCLVKGFYPSDIVVEWESSGQPENTYKTTPPVLDSDGSYFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSVSPGK"
  #mab = "mab1"
  check_input(hc_seq)
  hc_three <- to_aa_three_letter(hc_seq)
  hc_base <- c(hc_three, "water")
  hc_base_count_table <- count_molecules(hc_base)
  hc_to_elements <- get_seq_element_comp(hc_base)
  hc_elements_tibble <- hc_to_elements  %>%
    dplyr::left_join(hc_base_count_table, by = "molecule")

  # Cyclization check-----------------------------------------------------------
  hc_first_aa <- hc_three[1]
  if (hc_cyclized == TRUE & !(hc_first_aa %in% c("gln", "glu"))) {
    stop(paste0("HC can only be cyclized if it starts with Q (PyroQ) or E (PyroE)"),
         call. = FALSE)
  }

  # Clipping check--------------------------------------------------------------
  hc_last_aa <- hc_three[length(hc_three)]
  if (hc_clipped == TRUE & hc_last_aa != "lys") {
    stop(paste0("HC can only be clipped if ends with Lysine (K)"),
         call. = FALSE)
  }

  # Cyclization and clipping calculations of HC full reduced MW-----------------
  hc_first_aa <- hc_three[1]
  hc_last_aa <- hc_three[length(hc_three)]
  hc_cyclized = TRUE
  hc_clipped = TRUE

  if (hc_cyclized == TRUE & hc_clipped == TRUE) {
    if (hc_first_aa == "gln") {
      cli_cyc_tbl <- hc_elements_tibble %>%
        dplyr::bind_rows(
          tibble::tibble(
            element_composition %>% dplyr::filter(molecule %in% c("lys", "nh3")),
            n = -1L
          )
        )
    } else {
      cli_cyc_tbl <- hc_elements_tibble %>%
        dplyr::bind_rows(
          tibble::tibble(
            element_composition %>% dplyr::filter(molecule %in% c("lys", "water")),
            n = -1L
          )
        )
    }
  } else if (hc_cyclized == TRUE & hc_clipped == FALSE) {
    if (hc_first_aa == "gln") {
      cli_cyc_tbl <- hc_elements_tibble %>%
        dplyr::bind_rows(
          tibble::tibble(
            element_composition %>% dplyr::filter(molecule == "nh3"),
            n = -1L
          )
        )
    } else {
      cli_cyc_tbl <- hc_elements_tibble %>%
        dplyr::bind_rows(
          tibble::tibble(
            element_composition %>% dplyr::filter(molecule == "water"),
            n = -1L
          )
        )
    }
  } else if (hc_cyclized == FALSE & hc_clipped == TRUE) {
    cli_cyc_tbl <- hc_elements_tibble %>%
      dplyr::bind_rows(
        tibble::tibble(
          element_composition %>% dplyr::filter(molecule == "lys"),
          n = -1L
        )
      )
  } else if (hc_cyclized == FALSE & hc_clipped == FALSE) {
    cli_cyc_tbl <- hc_elements_tibble
  }

  hc_mab_full_reduced <- cli_cyc_tbl %>%
    dplyr::mutate(dplyr::across(2:(ncol(.)-1), ~.x*n)) %>%
    dplyr::select(-n) %>%
    dplyr::select(2:ncol(.)) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(.x, na.rm = TRUE))) %>%
    dplyr::mutate(mass = calculate_mass(.),
                  cyclyzation = ifelse(hc_cyclized == TRUE, "Yes", "No"),
                  clipping = ifelse(hc_clipped == TRUE, "Yes", "No")) %>%
    dplyr::mutate(antibody = mab, chain = "heavy", reduction = "full",
                  modifications = "No", glycans = "No", .before = 1)

  # Calculations for partially reduced HC---------------------------------------
  if (hc_cyclized == TRUE & hc_clipped == TRUE) {
    if (hc_first_aa == "gln") {
      partial_red_tbl <- hc_elements_tibble %>%
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
      partial_red_tbl <- hc_elements_tibble %>%
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
  } else if (hc_cyclized == TRUE & hc_clipped == FALSE) {
    if (hc_first_aa == "gln") {
      partial_red_tbl <- hc_elements_tibble %>%
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
      partial_red_tbl <- hc_elements_tibble %>%
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
  } else if (hc_cyclized == FALSE & hc_clipped == TRUE) {
    partial_red_tbl <- hc_elements_tibble %>%
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
  } else if (hc_cyclized == FALSE & hc_clipped == FALSE) {
    partial_red_tbl <- hc_elements_tibble %>%
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
    dplyr::mutate(mass = calculate_mass(.),
                  cyclyzation = ifelse(hc_cyclized == TRUE, "Yes", "No"),
                  clipping = ifelse(hc_clipped == TRUE, "Yes", "No")) %>%
    dplyr::mutate(antibody = mab, chain = "heavy", reduction = "partial",
                  modifications = "No", glycans = "No", .before = 1)

  # Add glycoforms to total reduced hc mass-------------------------------
  hc_glycan_tbl <- glycans %>%
    dplyr::filter(Show_HC == "Y") %>%
    dplyr::select(-c(Show_HC, Show_LC)) %>%
    tidyr::pivot_longer(-glycan_name,
                        names_to = "molecule",
                        values_to = "n")

  hc_glycan_individual <- split(hc_glycan_tbl , f = hc_glycan_tbl$glycan_name)

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
    dplyr::mutate(mass = calculate_mass(.),
                  cyclyzation = ifelse(hc_cyclized == TRUE, "Yes", "No"),
                  clipping = ifelse(hc_clipped == TRUE, "Yes", "No"),
                  antibody = mab, chain = "heavy", reduction = "full",
                  modifications = "No")

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
    dplyr::mutate(mass = calculate_mass(.),
                  cyclyzation = ifelse(hc_cyclized == TRUE, "Yes", "No"),
                  clipping = ifelse(hc_clipped == TRUE, "Yes", "No"),
                  antibody = mab, chain = "heavy", reduction = "partial",
                  modifications = "No")

  # Add user provided chem mod for HC

  # LC Mass calculation---------------------------------------------------------
  # lc_seq = "DIVLTQSPASLAVSLGQRATISCKASQSVDFDGDSFMNWYQQKPGQPPKLLIYTTSNLESGIPARFSASGSGTDFTLNIHPVEEEDTATYYCQQSNEDPYTFGGGTKLELKRAVAAPSVFIFPPSEDQVKSGTVSVVCLLNNFYPREASVKWKVDGVLKTGNSQESVTEQDSKDNTYSLSSTLTLSSTDYQSHNVYACEVTHQGLSSPVTKSFNRGEC"
  check_input(lc_seq)
  lc_three <- to_aa_three_letter(lc_seq)
  lc_base <- c(lc_three, "water")
  lc_base_count_table <- count_molecules(lc_base)
  lc_to_elements <- get_seq_element_comp(lc_base)
  lc_elements_tibble <- lc_to_elements  %>%
    dplyr::left_join(lc_base_count_table, by = "molecule")

  # Full reduction LC-----------------------------------------------------------
  lc_mab_full_reduced <- lc_elements_tibble %>%
    dplyr::mutate(dplyr::across(2:(ncol(.)-1), ~.x*n)) %>%
    dplyr::select(-n) %>%
    dplyr::select(2:ncol(.)) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(.x, na.rm = TRUE))) %>%
    dplyr::mutate(mass = calculate_mass(.),
                  cyclyzation = NA,
                  clipping = NA) %>%
    dplyr::mutate(antibody = mab, chain = "light", reduction = "full",
                  modifications = "No", glycans = "No", .before = 1)

  # Partial reduction LC--------------------------------------------------------
  partial_red_tbl_lc <- lc_elements_tibble %>%
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
    dplyr::mutate(mass = calculate_mass(.),
                  cyclyzation = NA,
                  clipping = NA) %>%
    dplyr::mutate(antibody = mab, chain = "light", reduction = "partial",
                  modifications = "No", glycans = "No", .before = 1)

  # Add glycoforms to total reduced lc mass-------------------------------------
  if (lc_glycosylation == T) {
    lc_glycan_tbl <- glycans %>%
      dplyr::filter(Show_LC == "Y") %>%
      dplyr::select(-c(Show_HC, Show_LC)) %>%
      tidyr::pivot_longer(-glycan_name,
                          names_to = "molecule",
                          values_to = "n")

    lc_glycan_individual <- split(lc_glycan_tbl , f = lc_glycan_tbl$glycan_name)

    lc_glycan_full <- lc_glycan_individual %>%
      purrr::map(dplyr::select, -1) %>%
      purrr::map(dplyr::left_join, element_composition, by = "molecule") %>%
      purrr::map(dplyr::relocate, n, .after = dplyr::last_col()) %>%
      purrr::map(dplyr::bind_rows, lc_elements_tibble) %>%
      purrr::map(dplyr::mutate_if, is.numeric, ~.x*n) %>%
      purrr::map(dplyr::select, -dplyr::last_col()) %>%
      purrr::map(dplyr::select, 2:dplyr::last_col()) %>%
      purrr::map(dplyr::summarise_all, ~ sum(.x, na.rm = TRUE))

    lc_glycan_full_calculation <- dplyr::bind_rows(lc_glycan_full) %>%
      dplyr::mutate(glycans = names(lc_glycan_full), .before = 1) %>%
      dplyr::mutate(mass = calculate_mass(.),
                    cyclyzation =NA,
                    clipping = NA,
                    antibody = mab, chain = "light", reduction = "full",
                    modifications = "No")

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
      dplyr::mutate(mass = calculate_mass(.),
                    cyclyzation =NA,
                    clipping = NA,
                    antibody = mab, chain = "light", reduction = "partial",
                    modifications = "No")
  } else {
    lc_glycan_full_calculation <- NULL
    lc_glycan_partial_calculation <- NULL
  }

  # Add user provided chem mod for LC

  return(dplyr::bind_rows(hc_mab_full_reduced,
                          hc_glycan_full_calculation,
                          hc_mab_partial_reduced,
                          hc_glycan_partial_calculation,
                          lc_mab_full_reduced,
                          lc_glycan_full_calculation,
                          lc_mab_partial_reduced,
                          lc_glycan_partial_calculation))
}
