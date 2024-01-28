#' Calculate charge and viscosity (in cP) of human IgG1 antibody at a pH of 5.5
#' Equation: η, cP(180mg/mL, 25°C) = 10^(0.15 + 1.26(0.60)∗ϕ − 0.043(0.047)∗q − 0.020(0.015)∗qsym)
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4284567/#eq1
#'
#' @param hc_fv A string containing variable heavy chain amino acids (one letter code)
#' @param lc_fv A string containing variable light chain amino acids (one letter code)
#'
#' @return A numeric value of predicted viscosity in cP
#' @export
#'
#' @examples
#' h_fv <- "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTR"
#' l_fv <- "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSAS"
#' visc <- mab_viscosity(h_fv, l_fv)
#' visc
#hc_fv = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVAS"
#lc_fv = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVERTVAA"
mab_viscosity <- function(hc_fv, lc_fv) {

  purrr::map(list(hc_fv, lc_fv), check_input)
  fv_charge <- Peptides::charge(stringr::str_c(hc_fv, lc_fv), pH = 5.5, pKscale = "Stryer")
  fv_csp <- Peptides::charge(hc_fv, pH = 5.5, pKscale = "Stryer") * Peptides::charge(lc_fv, pH = 5.5, pKscale = "Stryer")

  hi <- calc_hydrophobicity(stringr::str_c(hc_fv, lc_fv))

  return(tibble::tibble(fv_charge = fv_charge,
                        fvCSP = fv_csp,
                        hydrophobicity_index = hi,
                        viscosity = 10^(0.15 + (1.26 * 0.60 * hi) - (0.043 * 0.047 * fv_charge) - (0.020 * 0.015 * fv_csp))))

}
