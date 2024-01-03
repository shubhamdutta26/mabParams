# Create named vector for amino acids
aa_named <- c(A = "ala", R = "arg", N = "asn", D = "asp", C = "cys",
              Q = "glu", E = "gln", G = "gly", H = "his", I = "ile",
              L = "leu", K = "lys", M = "met", F = "phe", P = "pro",
              S = "ser", T = "thr", W = "trp", Y = "tyr", V = "val")

to_aa_three_letter <- function (one_letter_input) {
  one_letter_input %>%
    stringr::str_split(pattern = "", simplify = TRUE) %>%
    stringr::str_replace_all(aa_named)
}
