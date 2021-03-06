% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mkMutProfile.R
\name{mkMutProfileGeneCnv}
\alias{mkMutProfileGeneCnv}
\title{Add annotations to the gene cnv table}
\usage{
mkMutProfileGeneCnv(
  gene.cnv,
  arm.ploidies,
  min.cn.arm.ploidy.diff = 0.8,
  deep.del.max.max.copy.number = 0.3,
  trunc.max.min.copy.number = 0.3,
  loh.max.min.minor.allele.ploidy = 0.2,
  scoring = list(deep_deletion = c(5, 5), trunc = c(5, 5), loh = c(5, 0)),
  trunc.is.deep.del = T,
  verbose = T
)
}
\arguments{
\item{gene.cnv}{purple gene cnv dataframe}

\item{arm.ploidies}{arm ploidies dataframe}

\item{min.cn.arm.ploidy.diff}{Minimum difference between min_copy_number and arm_ploidy before a
gain is considered to have happened}

\item{deep.del.max.max.copy.number}{The max max_copy_number for a gene to be considered to 
be completely lost (deep deletion)}

\item{trunc.max.min.copy.number}{The max min_copy_number to for a gene to be considered to have
a truncation}

\item{loh.max.min.minor.allele.ploidy}{The maximum min_minor_allele_ploidy for a gene to be
considered to have loss of heterozygosity}

\item{scoring}{a list in the form list(deep_deletion=c(5,5),trunc=c(5,5),loh=c(5,0))}

\item{trunc.is.deep.del}{Consider truncations the same as deep deletions?}
}
\value{
The original input dataframe with annotation columns
}
\description{
Annotations deep deletions, truncations, LOH, and amplifications
}
