#' XXX
#'
#' @param hc_seq A string with hc amino acids
#' @param lc_seq A string with lc amino acids
#' @param hc_cyclized A boolean
#' @param hc_clipped A boolean
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
                        n_disulphides = 16L,
                        n_hc_disulphides = 4L,
                        n_lc_disulphides = 2L,
                        mab = "mab1") {
                       #lc_glycosylation = FALSE) {

  # Dataframes------------------------------------------------------------------
  # n = "n"
  # molecule = "molecule"
  element_composition <- tibble::tribble(
    ~molecule,~carbon,~hydrogen,~nitrogen,~oxygen,~sulphur,
    "carbon",1,0,0,0,0,
    "hydrogen",0,1,0,0,0,
    "nitrogen",0,0,1,0,0,
    "oxygen",0,0,0,1,0,
    "sulphur",0,0,0,0,1,
    "water",0,2,0,1,0,
    "nh3",0,3,1,0,0,
    "ala",3,5,1,1,0,
    "arg",6,12,4,1,0,
    "asn",4,6,2,2,0,
    "asp",4,5,1,3,0,
    "cys",3,5,1,1,1, # reduced with SH
    "glu",5,7,1,3,0,
    "gln",5,8,2,2,0,
    "gly",2,3,1,1,0,
    "his",6,7,3,1,0,
    "ile",6,11,1,1,0,
    "leu",6,11,1,1,0,
    "lys",6,12,2,1,0,
    "met",5,9,1,1,1,
    "phe",9,9,1,1,0,
    "pro",5,7,1,1,0,
    "ser",3,5,1,2,0,
    "thr",4,7,1,2,0,
    "trp",11,10,2,1,0,
    "tyr",9,9,1,2,0,
    "val",5,9,1,1,0,
    "HexNAc",8,13,1,5,0,
    "Hex",6,10,0,5,0,
    "dHex",6,10,0,4,0,
    "NeuAc",11,17,1,8,0,
    "NeuGc",11,17,1,9,0,
  )
  # HC Mass calculation---------------------------------------------------------
  # hc_seq = "QVTLKESGPGILQPSQTLSLTCSFSGFSLRTSGMGVGWIRQPSGKGLEWLAHIWWDDDKRYNPALKSRLTISKDTSSNQVFLKIASVDTADTATYYCAQINPAWFAYWGQGTLVTVSSASTKGPSVFPLAPSSRSTSESTAALGCLVKDYFPEPVTVSWNSGSLTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYVCNVNHKPSNTKVDKRVEIKTCGGGSKPPTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSQEDPDVKFNWYVNGAEVHHAQTKPRETQYNSTYRVVSVLTVTHQDWLNGKEYTCKVSNKALPAPIQKTISKDKGQPREPQVYTLPPSREELTKNQVSLTCLVKGFYPSDIVVEWESSGQPENTYKTTPPVLDSDGSYFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSVSPGK"
  # mab = "mab1"
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
  # hc_cyclized = TRUE
  # hc_clipped = TRUE

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



  # # Calculations for partially reduced HC
  # hc_mass_reduced_partial <- hc_mass_reduced_full - (n_hc_disulphides * (2 * 1.00794075))
  # adj_hc_reduced_partial <- cyc_cli(hc_mass_reduced_partial, hc_cyclized, hc_clipped, hc_first_aa)

  # Add user provided chem mod for HC
  # Add these glycoforms to total reduced hc mass

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

  # Add user provided chem mod for LC
  # Add these glycoforms to total reduced lc mass
  return(dplyr::bind_rows(hc_mab_full_reduced, hc_mab_partial_reduced, lc_mab_full_reduced, lc_mab_partial_reduced))
}
