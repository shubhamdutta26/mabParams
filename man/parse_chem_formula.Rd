% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_chem_formula.R
\name{parse_chem_formula}
\alias{parse_chem_formula}
\title{Converts chemical formula into a tibble with names and element counts. It does
not check if the letters match any chemical formula. For e.g. if X is present
in the formula. which is not in the periodic table it doesnot throw an error.
There cannot be any spaces or special characters in the formula.}
\usage{
parse_chem_formula(formula)
}
\arguments{
\item{formula}{A character vector og length 1 (e.g. C2H6O for ethanol; not C2H5OH)}
}
\value{
A tibble/ dataframe
}
\description{
Converts chemical formula into a tibble with names and element counts. It does
not check if the letters match any chemical formula. For e.g. if X is present
in the formula. which is not in the periodic table it doesnot throw an error.
There cannot be any spaces or special characters in the formula.
}
\examples{
formula <- "C6H12O6" # glucose
parse_chem_formula(formula)
}
