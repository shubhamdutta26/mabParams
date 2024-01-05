# Create named vector for amino acids
aa_one_letter <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                   "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
check_input <- function (one_letter_input) {
  if (typeof (one_letter_input) != "character" || length(one_letter_input) != 1) {
    stop("This needs to be a character vector of length 1", call. = FALSE)
  }
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

# Create named vector for amino acids
aa_named <- c(A = "ala", R = "arg", N = "asn", D = "asp", C = "cys",
              Q = "glu", E = "gln", G = "gly", H = "his", I = "ile",
              L = "leu", K = "lys", M = "met", F = "phe", P = "pro",
              S = "ser", T = "thr", W = "trp", Y = "tyr", V = "val")

to_aa_three_letter <- function (one_letter_input) {
  one_letter_input %>%
    stringr::str_split(pattern = "", simplify = TRUE) %>%
    toupper() %>%
    stringr::str_replace_all(aa_named)
}

# Count number of each molecule
count_molecules <- function (molecule) {
  table(molecule) %>%
    tibble::as_tibble()
}

# Elemental composition of each molecule
get_seq_element_comp <- function (seq) {

  element_composition <- tibble::tribble(
    ~molecule,~carbon,~hydrogen,~nitrogen,~oxygen,~sulphur,
    "carbon",1,0,0,0,0,
    "hydrogen",0,1,0,0,0,
    "nitrogen",0,0,1,0,0,
    "oxygen",0,0,0,1,0,
    "sulphur",0,0,0,0,1,
    "water",0,2,0,1,0,
    "ala",3,5,1,1,0,
    "arg",6,12,4,1,0,
    "asn",4,6,2,2,0,
    "asp",4,5,1,3,0,
    "cys",3,5,1,1,1,
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

  element_composition %>%
    dplyr::filter(element_composition$molecule %in% seq)
}

element_mass <- tibble::tribble(
  ~element,~symbol,~mass,
  "bromine","Br",79.904,
  "calcium","Ca",39.9625912,
  "carbon","C",12.0107359,
  "chlorine","Cl",35.4527,
  "fluorine","F",18.9984,
  "hydrogen","H",1.00794075,
  "iodine","I",126.90447,
  "lithium","Li",6.941,
  "nitrogen","N",14.00670321,
  "oxygen","O",15.99940492,
  "phosphorus","P",30.97376,
  "potassium","K",39.0983,
  "selenium","Se",78.96,
  "sodium","Na",22.98977,
  "sulphur","S",32.06478741,
)

element_mass_list <- list(
  bromine=79.904,
  calcium=39.9625912,
  carbon=12.0107359,
  chlorine=35.4527,
  fluorine=18.9984,
  hydrogen=1.00794075,
  iodine=126.90447,
  lithium=6.941,
  nitrogen=14.00670321,
  oxygen=15.99940492,
  phosphorus=30.97376,
  potassium=39.0983,
  selenium=78.96,
  sodium=22.98977,
  sulphur=32.06478741
)

calculate_mass <- function (data, n, counts, mass, each_mass) {
  data %>%
    dplyr::mutate(dplyr::across(2:(ncol(data)-1), ~.x*n)) %>%
    dplyr::select(-n) %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(names(element_mass_list)),
                                ~.x * element_mass_list[[dplyr::cur_column()]])) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sum_each_mol = sum(dplyr::c_across(dplyr::where(is.numeric)))) %>%
    dplyr::pull("sum_each_mol") %>%
    sum()
}