% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TFBS.findR.R
\name{p.value.calculation}
\alias{p.value.calculation}
\title{p.value.calculation}
\usage{
p.value.calculation(data, pwms, background = c(A = 0.25, C = 0.25, G =
  0.25, T = 0.25))
}
\arguments{
\item{data}{is a Grangesobject from \code{update.data.R}}

\item{pwms}{is a simplelist of motifs}

\item{background}{is the background probability for each nucleotide.}
}
\value{
Grangesobject with added columns:
\item{Pvalue.ref.score}{Pvalue of transcription factor binding score
for each position and motif to Ref.sequence}
\item{Pvalue.alt.score}{Pvalue of transcription factor binding score
for each position and motif to Alt.sequence}
}
\description{
Touzet et al., 2007, created a method for calculating
pvalues for PWMs
https://almob.biomedcentral.com/articles/10.1186/1748-7188-2-15
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0480-5
}
