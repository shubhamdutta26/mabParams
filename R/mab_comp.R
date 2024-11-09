#' mab_comp returns a tibble with the amino acid compositions of intact
#' monoclonal antibody (cys residues are reduced). Charges on each amino
#' acid is assigned based on pI of the antibody
#'
#' @param hc_seq A string containing heavy chain amino acids (one letter code)
#' @param lc_seq A string containing light chain amino acids (one letter code)
#'
#' @return A tibble with amino acid compositions
#' @export
#'
#' @examples
#' heavy <- "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYA"
#' light <- "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSR"
#' comp_tbl <- mab_comp(heavy, light)
#' head(comp_tbl)
mab_comp <- function (hc_seq, lc_seq) {
  purrr::map(list(hc_seq, lc_seq), check_input)
  hc <- to_aa_three_letter(hc_seq)
  lc <- to_aa_three_letter(lc_seq)
  hc_count <- count_molecules(hc)%>%
    dplyr::mutate(chain = "heavy", .before = 1) %>%
    dplyr::rename(amino_acid = "molecule") %>%
    dplyr::mutate(n = n*2,
                  aa_percent = round(n/sum(n), 2),
                  aa_type = dplyr::case_match(
                    amino_acid,
                    c("ala","cys","gly","ser","thr") ~ "tiny",
                    c("ala","cys","asp","gly","asn","pro","ser","thr","val") ~ "small",
                    c("ala","ile","leu","val") ~ "aliphatic",
                    c("phe","his","trp","tyr") ~ "aromatic",
                    c("ala","cys","phe","gly","ile","leu","met","pro","val","trp","tyr") ~ "non-polar",
                    c("asp","glu","his","lys","asn","gln","arg","ser","thr") ~ "polar",
                    c("asp","glu","his","lys","arg") ~ "charged",
                    c("his","lys","arg") ~ "basic",
                    c("asp","glu") ~ "acidic"
                  ))

  lc_count <- count_molecules(lc)%>%
    dplyr::mutate(chain = "light", .before = 1) %>%
    dplyr::rename(amino_acid = "molecule") %>%
    dplyr::mutate(n = n*2,
                  aa_percent = round(n/sum(n), 2),
                  aa_type = dplyr::case_match(
                    amino_acid,
                    c("ala","cys","gly","ser","thr") ~ "tiny",
                    c("ala","cys","asp","gly","asn","pro","ser","thr","val") ~ "small",
                    c("ala","ile","leu","val") ~ "aliphatic",
                    c("phe","his","trp","tyr") ~ "aromatic",
                    c("ala","cys","phe","gly","ile","leu","met","pro","val","trp","tyr") ~ "non-polar",
                    c("asp","glu","his","lys","asn","gln","arg","ser","thr") ~ "polar",
                    c("asp","glu","his","lys","arg") ~ "charged",
                    c("his","lys","arg") ~ "basic",
                    c("asp","glu") ~ "acidic"
                  ))

  all_count <- count_molecules(c(hc, lc))%>%
    dplyr::mutate(chain = "both", .before = 1) %>%
    dplyr::rename(amino_acid = "molecule") %>%
    dplyr::mutate(n = n*2,
                  aa_percent = round(n/sum(n), 2),
                  aa_type = dplyr::case_match(
                    amino_acid,
                    c("ala","cys","gly","ser","thr") ~ "tiny",
                    c("ala","cys","asp","gly","asn","pro","ser","thr","val") ~ "small",
                    c("ala","ile","leu","val") ~ "aliphatic",
                    c("phe","his","trp","tyr") ~ "aromatic",
                    c("ala","cys","phe","gly","ile","leu","met","pro","val","trp","tyr") ~ "non-polar",
                    c("asp","glu","his","lys","asn","gln","arg","ser","thr") ~ "polar",
                    c("asp","glu","his","lys","arg") ~ "charged",
                    c("his","lys","arg") ~ "basic",
                    c("asp","glu") ~ "acidic"
                  ))

  return(dplyr::bind_rows(all_count, hc_count, lc_count))
}
