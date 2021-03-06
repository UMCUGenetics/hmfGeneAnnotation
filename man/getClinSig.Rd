% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preProcessHmfPlOutput.R
\name{getClinSig}
\alias{getClinSig}
\title{Get clinical significance by querying a tabix indexed bed file}
\usage{
getClinSig(df, db.path, verbose = T)
}
\arguments{
\item{df}{A bed file like dataframe containing the columns: chrom, pos, ref, alt}

\item{db.path}{Path to a tabix indexed bed file with the clinical significance annotations (e.g.
from ClinVar). This bed file should contain the columns (in this order): chrom, pos, ref, alt, 
clinsig, id. However, it should not have colnames.}

\item{verbose}{Show progress messages?}
}
\value{
A character vector of the clinical significance annotations
}
\description{
Get clinical significance by querying a tabix indexed bed file
}
