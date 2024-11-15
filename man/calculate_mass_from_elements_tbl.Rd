% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{calculate_mass_from_elements_tbl}
\alias{calculate_mass_from_elements_tbl}
\title{Calculate Molecular Mass from Element Composition}
\usage{
calculate_mass_from_elements_tbl(df)
}
\arguments{
\item{df}{A data frame where column names match element names (e.g., "carbon", "hydrogen")
and values represent the number of atoms of each element}
}
\value{
A numeric vector containing the calculated molecular mass for each row
}
\description{
Takes a data frame containing element quantities and calculates the total molecular mass
for each row by multiplying element quantities with their respective atomic masses
and summing the results.
}
\details{
The function uses a predefined \code{element_mass_list} containing atomic masses for elements.
Only columns matching element names in this list are used in calculations.
Non-matching columns are ignored. Missing values (NA) are treated as 0 in the sum.
}
\examples{
\dontrun{
compounds <- data.frame(
  carbon = c(1, 2),
  hydrogen = c(4, 6),
  oxygen = c(1, 1)
)
calculate_mass_from_elements_tbl(compounds)
}

}
\keyword{internal}
