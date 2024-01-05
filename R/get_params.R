#' XXX
#'
#' @param hc_seq A string with hc amino acids
#' @param lc_seq A string with lc amino acids
#' @return vector
#' @export
#'
#' @examples
get_params_monoclonal <- function (hc_seq,
                                   lc_seq) {
                                   #hc_cyclized = TRUE,
                                   #hc_clipped = TRUE,
                                   #n_disulphides = 16,
                                   #unred_hc_disulphides = 4,
                                   #unred_lc_disulphides = 2,
                                   #lc_glycosylation = FALSE) {
  # For dplyr-------------------------------------------------------------------
  n = "n"
  counts = "counts"
  mass = "mass"
  each_mass = "each_mass"
  symbol = "symbol"

  # HC Mass calculation---------------------------------------------------------
  hc_seq = "QVTLKESGPGILQPSQTLSLTCSFSGFSLRTSGMGVGWIRQPSGKGLEWLAHIWWDDDKRYNPALKSRLTISKDTSSNQVFLKIASVDTADTATYYCAQINPAWFAYWGQGTLVTVSSASTKGPSVFPLAPSSRSTSESTAALGCLVKDYFPEPVTVSWNSGSLTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYVCNVNHKPSNTKVDKRVEIKTCGGGSKPPTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSQEDPDVKFNWYVNGAEVHHAQTKPRETQYNSTYRVVSVLTVTHQDWLNGKEYTCKVSNKALPAPIQKTISKDKGQPREPQVYTLPPSREELTKNQVSLTCLVKGFYPSDIVVEWESSGQPENTYKTTPPVLDSDGSYFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSVSPGK"
  check_input(hc_seq)
  aa_three_hc <- to_aa_three_letter(hc_seq)
  count_table_hc <- count_molecules(aa_three_hc)
  element_table_hc <- get_seq_element_comp(aa_three_hc)
  joined_hc <- element_table_hc %>%
    dplyr::left_join(count_table_hc, by = "molecule")
  hc_mass <- calculate_mass(joined_hc, n, counts, mass, each_mass)
  # hc_mass
  # Cyclization check

  # nh3_mass <- element_mass %>%
  #   dplyr::filter(symbol == "N") %>%
  #   dplyr::pull(mass) +
  #   (element_mass %>%
  #      dplyr::filter(symbol == "H") %>%
  #      dplyr::pull(mass))*3
  #
  # h2o_mass <- element_mass %>%
  #   dplyr::filter(symbol == "H") %>%
  #   dplyr::pull(mass) +
  #   (element_mass %>%
  #      dplyr::filter(symbol == "O") %>%
  #      dplyr::pull(mass))*2
  #
  # if (hc_cyclized == TRUE) {
  #   if (stringr::str_starts(hc_seq, "Q") == TRUE) {
  #     hc_mass = hc_mass - nh3_mass
  #   } else if (stringr::str_starts(hc_seq, "E") == TRUE) {
  #     hc_mass = hc_mass - h2o_mass
  #   } else {
  #     stop(paste0("HC can only be cyclized if it starts with Q (PyroQ) or E (PyroE)"),
  #          call. = FALSE)
  #   }
  # }
  #
  # # Lysine clipping check
  # K_mass = (6*12.0107359)+(12*1.00794075)+(2*14.00670321)+(1*15.99940492)
  #
  # if (hc_clipped == TRUE) {
  #   if (stringr::str_ends(hc_seq, "K") == TRUE) {
  #     hc_mass = hc_mass - K_mass
  #   }
  # }
  #
  # # For fully reduced mass 11
  # full_reduced_hc_mass <- hc_mass + (11 * 1.00794075)
  # partial_reduced_hc_mass <- full_reduced_hc_mass - (unred_hc_disulphides*2*1.00794075)
  #
  # # Add user provided chem mod for HC
  # # Add these glycoforms to total reduced hc mass
  #
  # # LC Mass calculation---------------------------------------------------------
  lc_seq = "DIVLTQSPASLAVSLGQRATISCKASQSVDFDGDSFMNWYQQKPGQPPKLLIYTTSNLESGIPARFSASGSGTDFTLNIHPVEEEDTATYYCQQSNEDPYTFGGGTKLELKRAVAAPSVFIFPPSEDQVKSGTVSVVCLLNNFYPREASVKWKVDGVLKTGNSQESVTEQDSKDNTYSLSSTLTLSSTDYQSHNVYACEVTHQGLSSPVTKSFNRGEC"
  check_input(lc_seq)
  aa_three_lc <- to_aa_three_letter(lc_seq)
  count_table_lc <- count_molecules(aa_three_lc)
  element_table_lc <- get_seq_element_comp(aa_three_lc)
  joined_lc <- element_table_lc %>%
    dplyr::left_join(count_table_lc, by = "molecule")
  lc_mass <- calculate_mass(joined_lc, n, counts, mass, each_mass)
  #lc_mass

  (hc_mass + lc_mass)*2

}
