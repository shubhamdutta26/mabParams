# Create named vector for amino acids
aa_one_letter <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                   "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

check_input <- function (one_letter_input) {
  split_aa <- one_letter_input %>%
    stringr::str_split(pattern = "", simplify = TRUE) %>%
    as.character() %>%
    toupper()
  indexes = !(split_aa %in% aa_one_letter)
  if (sum(indexes) >= 1) {
    x <- stringr::str_flatten_comma(split_aa[indexes == TRUE])
    stop(paste0("Input contains :", x, "; Only one letter amino acid codes should be in the input without spaces."), call. = FALSE)
  }
}
