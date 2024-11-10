#' Calculate Amino Acid Composition of Monoclonal Antibody Sequences
#'
#' @description
#' Analyzes the amino acid composition of heavy and light chain sequences of a monoclonal
#' antibody (MAb). For individual chains, counts represent single chain composition.
#' For both chains combined, counts are doubled to represent the complete antibody
#' (2 heavy + 2 light chains).
#'
#' @param hc_seq Character string of single-letter amino acid codes for heavy chain
#' @param lc_seq Character string of single-letter amino acid codes for light chain
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{chain}{Character indicating "both", "heavy", or "light" chain}
#'   \item{amino_acid}{Three-letter amino acid code}
#'   \item{n}{Count of amino acid occurrences (doubled only for "both" to represent complete antibody)}
#'   \item{aa_percent}{Percentage of each amino acid (rounded to 2 decimal places)}
#'   \item{polarity}{Amino acid polarity classification}
#'   \item{net_charge_7.4}{Net charge at pH 7.4}
#' }
#'
#' @importFrom dplyr mutate rename bind_rows left_join select
#' @importFrom tibble tribble
#'
#' @examples
#' hc <- "QVQLVQSGAEVKKPGASVKVSCKASGYTFT"
#' lc <- "DIQMTQSPSSLSASVGDRVTITCRASQ"
#' mab_comp(hc, lc)
#'
#' @export
mab_comp <- function(hc_seq, lc_seq) {
  # Input validation
  check_input(hc_seq)
  check_input(lc_seq)

  # Convert to three-letter codes
  hc <- to_aa_three_letter(hc_seq)
  lc <- to_aa_three_letter(lc_seq)

  # Helper function to process counts
  process_counts <- function(seq, chain_name, double_counts = FALSE) {
    counts <- count_molecules(seq) %>%
      dplyr::mutate(
        chain = chain_name,
        .before = 1
      ) %>%
      dplyr::rename(amino_acid = "molecule")

    # Only double counts for complete antibody (both chains)
    if(double_counts) {
      counts <- counts %>% dplyr::mutate(n = n * 2)
    }

    counts %>%
      dplyr::mutate(
        aa_percent = round(n/sum(n), 3)
      ) %>%
      # Join with aa_type table for properties
      dplyr::left_join(
        aa_type %>%
          dplyr::rename(amino_acid = aa),
        by = "amino_acid"
      )
  }

  # Process individual chains without doubling
  hc_count <- process_counts(hc, "heavy", double_counts = FALSE)
  lc_count <- process_counts(lc, "light", double_counts = FALSE)

  # Process both chains together with doubling
  all_count <- process_counts(c(hc, lc), "both", double_counts = TRUE)

  # Combine results
  dplyr::bind_rows(all_count, hc_count, lc_count)
}
