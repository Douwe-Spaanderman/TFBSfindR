% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TFBS.findR.R
\name{to.PWM}
\alias{to.PWM}
\title{to.PWM}
\usage{
to.PWM(x, type = "PFM", background = c(A = 0.25, C = 0.25, G = 0.25, T
  = 0.25), pseudocount = 0.8)
}
\arguments{
\item{x}{is a GrangesObject from \code{read.input.file}}

\item{type}{is the type of input file, which can be "PFM", "PCM",
"PPM" and "PWM". If the latter is passed, than the PWM is returned
without changes}

\item{background}{is the background probability for each nucleotide.}

\item{pseudocount}{value to change the outcome model, which removes
possibility for -inf values.}
}
\value{
transformed matrix of x based on the method selected
}
\description{
Function to change position frequency matrix/position count matrix or
position probability matrix to a position weight matrix.
}
\details{
Several methods have been discussed to calculate PWM,
Papatsenko et al., 2003, suggested the method we use now, with the use of
a pseudocount and sum log. As suggested by papatsenko, pseudocount is
calculated by log(sequenceCount), in which sequenceCount is the amount
of times PFM is sequenced. If however this information isn't available
-or the other method is prefered- than the pseudocount can be set
beforehand. Nakai et al., 2008, discussed than smaller pseudocounts
preserve data model better than larger pseudocount. The lowest value
they analysed was 0.01 and therefore this is also our default.
Note that if you want to calculate the pseudocount based on the log of
sequence reads, this param should be set to log.of.reads.
}
