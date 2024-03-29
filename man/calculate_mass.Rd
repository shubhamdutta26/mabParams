% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_mass.R
\name{calculate_mass}
\alias{calculate_mass}
\title{Calculate mass of peptides (reduced) in daltons}
\usage{
calculate_mass(seq)
}
\arguments{
\item{seq}{One peptide sequence without spaces}
}
\value{
The numeric mass in daltons
}
\description{
Calculate mass of peptides (reduced) in daltons
}
\examples{
sequence <- "QVTLKESGPGILQPSQTLSLTCSFSGFSLRTSGMGV"
mass <- calculate_mass(sequence)
print(mass)
}
