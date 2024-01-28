# Create named vector for amino acids
# aa_one_letter <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
#                    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

# Check validity of HC/LC input
check_input <- function (one_letter_input) {

  if (typeof (one_letter_input) != "character" || length(one_letter_input) != 1) {
    stop("HC/ LC sequence input needs to be a character vector of length 1", call. = FALSE)
  }

  split_aa <- toupper(as.character(stringr::str_split(one_letter_input, pattern = "", simplify = TRUE)))

  indexes = !(split_aa %in% aa_one_letter)

  if (sum(indexes) >= 1) {
    x <- stringr::str_flatten_comma(split_aa[indexes == TRUE])
    stop(paste0("Input contains:", x, "; Only one letter amino acid codes should be in the input without spaces."), call. = FALSE)
  }
}

check_boolean_args <- function(argument) {
  if (typeof(argument) != "logical"){
    stop(paste0("Cyclization, clipping, and glycosylation arguments must be logical."))
  }
}

check_chem_mod_args <- function(hc_chem_mod = NA, lc_chem_mod = NA) {
  if (!is.na(hc_chem_mod)) {
    if (typeof(hc_chem_mod) != "character" || length(hc_chem_mod) > 1) {
      stop("The chemical modifications needs to be a character vector of length 1.", call. = FALSE)
    }
    if (grepl("^[0-9]$", substr(hc_chem_mod, 1, 1)) == TRUE || grepl("[^a-zA-Z0-9]", hc_chem_mod) == TRUE ) {
      stop("The chemical modifications cannot begin with a number or contain special characters.", call. = FALSE)
    }
  }
  if (!is.na(lc_chem_mod)) {
    if (typeof(lc_chem_mod) != "character" || length(lc_chem_mod) > 1) {
      stop("The chemical modifications needs to be a character vector of length 1.", call. = FALSE)
    }
    if (grepl("^[0-9]$", substr(lc_chem_mod, 1, 1)) == TRUE || grepl("[^a-zA-Z0-9]", lc_chem_mod) == TRUE ) {
      stop("The chemical modifications cannot begin with a number or contain special characters.", call. = FALSE)
    }
  }
}

# Create named vector for amino acids
# aa_named <- c(A = "ala", R = "arg", N = "asn", D = "asp", C = "cys",
#               Q = "gln", E = "glu", G = "gly", H = "his", I = "ile",
#               L = "leu", K = "lys", M = "met", F = "phe", P = "pro",
#               S = "ser", T = "thr", W = "trp", Y = "tyr", V = "val")

# Convert HC and LC to three letter codes
to_aa_three_letter <- function (one_letter_input) {
  toupper(stringr::str_split(one_letter_input, pattern = "", simplify = TRUE)) %>%
    stringr::str_replace_all(aa_named)
}

# Count number of each molecule
count_molecules <- function (molecule) {
  tibble::as_tibble(table(molecule))
}

# Elemental composition of each molecule
get_seq_element_comp <- function (seq) {

  # element_composition <- tibble::tribble(
  #   ~molecule,~carbon,~hydrogen,~nitrogen,~oxygen,~sulphur,
  #   "carbon",1,0,0,0,0,
  #   "hydrogen",0,1,0,0,0,
  #   "nitrogen",0,0,1,0,0,
  #   "oxygen",0,0,0,1,0,
  #   "sulphur",0,0,0,0,1,
  #   "water",0,2,0,1,0,
  #   "nh3",0,3,1,0,0,
  #   "ala",3,5,1,1,0,
  #   "arg",6,12,4,1,0,
  #   "asn",4,6,2,2,0,
  #   "asp",4,5,1,3,0,
  #   "cys",3,5,1,1,1, # reduced with SH
  #   "glu",5,7,1,3,0,
  #   "gln",5,8,2,2,0,
  #   "gly",2,3,1,1,0,
  #   "his",6,7,3,1,0,
  #   "ile",6,11,1,1,0,
  #   "leu",6,11,1,1,0,
  #   "lys",6,12,2,1,0,
  #   "met",5,9,1,1,1,
  #   "phe",9,9,1,1,0,
  #   "pro",5,7,1,1,0,
  #   "ser",3,5,1,2,0,
  #   "thr",4,7,1,2,0,
  #   "trp",11,10,2,1,0,
  #   "tyr",9,9,1,2,0,
  #   "val",5,9,1,1,0,
  #   "HexNAc",8,13,1,5,0,
  #   "Hex",6,10,0,5,0,
  #   "dHex",6,10,0,4,0,
  #   "NeuAc",11,17,1,8,0,
  #   "NeuGc",11,17,1,9,0,
  # )


  element_composition %>%
    dplyr::filter(element_composition$molecule %in% seq)
}

# element_mass_list <- list(
#   bromine=79.904,
#   calcium=39.9625912,
#   carbon=12.0107359,
#   chlorine=35.4527,
#   fluorine=18.9984,
#   hydrogen=1.00794075,
#   iodine=126.90447,
#   lithium=6.941,
#   nitrogen=14.00670321,
#   oxygen=15.99940492,
#   phosphorus=30.97376,
#   potassium=39.0983,
#   selenium=78.96,
#   sodium=22.98977,
#   sulphur=32.06478741
# )

calculate_mass_from_elements_tbl <- function (df) {
  df %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(names(element_mass_list)),
    ~.x * element_mass_list[[dplyr::cur_column()]])) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(mass_full = sum(dplyr::c_across(dplyr::where(is.numeric)))) %>%
    dplyr::pull("mass_full")
}

combine_and_sum <- function(df1, df2, char_col, num_cols, suffix_1, suffix_2) {
  # Create all combinations of the 'char' column
  comb <- expand.grid(df1[[char_col]], df2[[char_col]])

  # Name the columns
  names(comb) <- c(paste0(char_col, suffix_1), paste0(char_col, suffix_2))

  # Loop over the numeric columns
  for (num_col in num_cols) {
    comb[[paste0(num_col, suffix_1)]] <- df1[[num_col]][match(comb[[paste0(char_col, suffix_1)]], df1[[char_col]])]
    comb[[paste0(num_col, suffix_2)]] <- df2[[num_col]][match(comb[[paste0(char_col, suffix_2)]], df2[[char_col]])]
    comb[[paste0("sum_", num_col)]] <- comb[[paste0(num_col, suffix_1)]] + comb[[paste0(num_col, suffix_2)]]
  }

  return(comb)
}

calc_hydrophobicity <- function(seq) {
  # https://www.sciencedirect.com/science/article/abs/pii/0022283684903097?fr=RR-2&ref=pdf_download&rr=84c215f36a7b4d19
  # aa_hydrophobic <- readr::read_csv(file = "inst/extdata/aa_hydrophobicity.csv",
  #                                   col_types = "ccnnc")
  tbl <- count_molecules(to_aa_three_letter(seq)) %>%
    dplyr::rowwise() %>%
    dplyr::full_join(aa_hydrophobic, by = dplyr::join_by(molecule == aa_three)) %>%
    dplyr::mutate(product = n * norm_consensus)
  hi <- -(sum(subset(tbl, hydrophobic == "yes")$product, na.rm = TRUE)/sum(subset(tbl, hydrophobic == "no")$product, na.rm = TRUE))
  return(hi)
}
