#' XXX
#'
#' @param hc_seq A string with hc amino acids
#' @param lc_seq A string with lc amino acids
#' @param hc_cyclized A boolean
#' @param hc_clipped A boolean
#' @param n_disulphides Number of total disulphides
#' @param n_hc_disulphides Number of total HC disulphides
#' @param n_lc_disulphides Number of total LC disulphides
#' @return vector
#' @export
#'
#' @examples
mab_params <- function (hc_seq,
                        lc_seq,
                        hc_cyclized,
                        hc_clipped,
                        n_disulphides = 16,
                        n_hc_disulphides = 4,
                        n_lc_disulphides = 2) {
                       #lc_glycosylation = FALSE) {

  # HC Mass calculation---------------------------------------------------------
  # hc_seq = "QVTLKESGPGILQPSQTLSLTCSFSGFSLRTSGMGVGWIRQPSGKGLEWLAHIWWDDDKRYNPALKSRLTISKDTSSNQVFLKIASVDTADTATYYCAQINPAWFAYWGQGTLVTVSSASTKGPSVFPLAPSSRSTSESTAALGCLVKDYFPEPVTVSWNSGSLTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYVCNVNHKPSNTKVDKRVEIKTCGGGSKPPTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSQEDPDVKFNWYVNGAEVHHAQTKPRETQYNSTYRVVSVLTVTHQDWLNGKEYTCKVSNKALPAPIQKTISKDKGQPREPQVYTLPPSREELTKNQVSLTCLVKGFYPSDIVVEWESSGQPENTYKTTPPVLDSDGSYFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSVSPGK"
  check_input(hc_seq)
  aa_three_hc <- to_aa_three_letter(hc_seq)
  count_table_hc <- count_molecules(aa_three_hc)
  element_table_hc <- get_seq_element_comp(aa_three_hc)
  joined_hc <- element_table_hc %>%
    dplyr::left_join(count_table_hc, by = "molecule")
  hc_mass <- calculate_mass(joined_hc, n = "n")
  hc_mass_reduced_full <- hc_mass + 18.01528 # Add water

  # Cyclization check-----------------------------------------------------------
  hc_first_aa <- aa_three_hc[1]
  if (hc_cyclized == TRUE & !(hc_first_aa %in% c("gln", "glu"))) {
    stop(paste0("HC can only be cyclized if it starts with Q (PyroQ) or E (PyroE)"),
         call. = FALSE)
  }

  # Clipping check--------------------------------------------------------------
  hc_last_aa <- aa_three_hc[length(aa_three_hc)]
  if (hc_clipped == TRUE & hc_last_aa != "lys") {
    stop(paste0("HC can only be clipped if ends with Lysine (K)"),
         call. = FALSE)
  }

  # Cyclization and clipping check and calculations of HC full reduced MW
  adj_hc_reduced_full <- cyc_cli(hc_mass_reduced_full, hc_cyclized, hc_clipped, hc_first_aa)

  # Calculations for partially reduced HC
  hc_mass_reduced_partial <- hc_mass_reduced_full - (n_hc_disulphides * (2 * 1.00794075))
  adj_hc_reduced_partial <- cyc_cli(hc_mass_reduced_partial, hc_cyclized, hc_clipped, hc_first_aa)

  # Add user provided chem mod for HC
  # Add these glycoforms to total reduced hc mass

  # LC Mass calculation---------------------------------------------------------
  # lc_seq = "DIVLTQSPASLAVSLGQRATISCKASQSVDFDGDSFMNWYQQKPGQPPKLLIYTTSNLESGIPARFSASGSGTDFTLNIHPVEEEDTATYYCQQSNEDPYTFGGGTKLELKRAVAAPSVFIFPPSEDQVKSGTVSVVCLLNNFYPREASVKWKVDGVLKTGNSQESVTEQDSKDNTYSLSSTLTLSSTDYQSHNVYACEVTHQGLSSPVTKSFNRGEC"
  check_input(lc_seq)
  aa_three_lc <- to_aa_three_letter(lc_seq)
  count_table_lc <- count_molecules(aa_three_lc)
  element_table_lc <- get_seq_element_comp(aa_three_lc)
  joined_lc <- element_table_lc %>%
    dplyr::left_join(count_table_lc, by = "molecule")
  lc_mass <- calculate_mass(joined_lc, n = "n")
  lc_mass_reduced_full <- lc_mass + 18.01528 # Add water

  # Calculations for partially reduced LC
  lc_mass_reduced_partial <- lc_mass_reduced_full - (n_lc_disulphides * (2 * 1.00794075))

  # Add user provided chem mod for LC
  # Add these glycoforms to total reduced lc mass

  # Intact Mab calculation------------------------------------------------------
  mab_mass_intact <- ((adj_hc_reduced_full + lc_mass_reduced_full) * 2) - (n_disulphides * 2 * 1.00794075)

  # Return base masses with no glycosylation or modifications
  return(list(mab_mass_intact, adj_hc_reduced_full, lc_mass_reduced_full,
              adj_hc_reduced_partial, lc_mass_reduced_partial))

  # Add user provided chem mod for LC
  # Add these glycoforms to total reduced lc mass
}
