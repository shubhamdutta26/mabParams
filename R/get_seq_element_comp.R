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
