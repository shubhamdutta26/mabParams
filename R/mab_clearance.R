#' Predicts clearance of human IgG1 based on variable fragments and CDRs of
#' heavy and light chains
#'
#' @param hc_fv A string containing heavy variable fragment amino acids (one letter code)
#' @param lc_fv A string containing light variable fragment amino acids (one letter code)
#' @param cdr_h3 A string containing heavy CDR3 amino acids (one letter code)
#' @param cdr_l1 A string containing heavy CDL1 amino acids (one letter code)
#' @param cdr_l3 A string containing heavy CDL1 amino acids (one letter code)
#'
#' @return A tibble with clearance parameters
#' @export
#'
#' @examples
#' h_fv <- "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTR"
#' l_fv <- "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSAS"
#' cdr_h3 = "WGGDGFYAMDY"
#' cdr_l1 = "RASQDVNTAVA"
#' cdr_l3 = "SASFLYS"
#' clearance_df <- mab_clearance(h_fv, l_fv, cdr_h3, cdr_l1, cdr_l3)
#' clearance_df
mab_clearance <- function(hc_fv, lc_fv, cdr_h3, cdr_l1, cdr_l3) {
  #hc_fv = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVAS"
  #lc_fv = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVERTVAA"
  #cdr_h3 = "WGGDGFYAMDY"
  #cdr_l1 = "RASQDVNTAVA"
  #cdr_l3 = "SASFLYS"
  purrr::map(list(hc_fv, lc_fv, cdr_h3, cdr_l1, cdr_l3), check_input)

  fv_charge <- Peptides::charge(stringr::str_c(hc_fv, lc_fv), pH = 5.5, pKscale = "Stryer")

  sum_hi <- calc_hydrophobicity(cdr_h3) + calc_hydrophobicity(cdr_l1) + calc_hydrophobicity(cdr_l3)

  final_tbl <- tibble::tibble(
    fv_charge = fv_charge,
    sum_hydrophobicity_index = sum_hi) %>%
    dplyr::mutate(
      clearance =  dplyr::case_when(
        (fv_charge >= 0 & fv_charge <= 6.2 & sum_hydrophobicity_index <= 4) ~ "Normal (< 10ml/kg/day)",
        (fv_charge < 0 & fv_charge > 6.2 & sum_hydrophobicity_index > 4) ~ "Fast (>= 10ml/kg/day)"
      )
    )

  return(final_tbl)
}
