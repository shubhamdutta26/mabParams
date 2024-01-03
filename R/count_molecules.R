# Count number of each molecule
count_molecules <- function (molecule) {
  table(molecule) %>%
    tibble::as_tibble()
}
