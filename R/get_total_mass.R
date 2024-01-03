get_total_mass <- function (df) {

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

  df %>%
    dplyr::mutate(dplyr::across(2:(ncol(df)-1), ~.x*n)) %>%
    dplyr::select(-n) %>%
    tidyr::pivot_longer(-1, names_to = "element", values_to = "counts") %>%
    dplyr::left_join(element_mass, by = "element") %>%
    dplyr::mutate(each_mass = counts * mass) %>%
    dplyr::pull(each_mass) %>%
    sum()
}
