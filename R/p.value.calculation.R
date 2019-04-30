#' p.value.calculation
#'
#' Touzet et al., 2007, created a method for calculating
#' pvalues for PWMs
#'
#' @param data is a Grangesobject from \code{update.data.R}
#' @param motiflist is a simplelist of motifs
#' @param background is the postion frequency matrix with each
#' background score of a nucleotide
#' @param motif.type is the type of motifs passed. These can be
#' "PWM", "PFM", "PCM" or "PPM".
#' @param pseudocount value to change the outcome model, which removes
#' possibility for -inf values.
#'
#' @return Grangesobject with added columns:
#' \item{Pvalue.ref.score}{Pvalue of transcription factor binding score
#' for each position and motif to Ref.sequence}
#' \item{Pvalue.alt.score}{Pvalue of transcription factor binding score
#' for each position and motif to Alt.sequence}
#'
#' @import TFMPvalue
#' @import S4Vectors
#' @export
p.value.calculation <- function(data, motiflist, background=c(A=0.25, C=0.25, G=0.25, T=0.25), motif.type="PFM", pseudocount=0.01){
  #Add names to pwms
  for (i in 1:length(motiflist)){
    names(motiflist[[i]]) <- names(motiflist[i])
    names(motiflist[[i]])[2] <- mcols(motiflist)$sequenceCount[i]
  }
  #Changing motifs to PWMs
  pwms <- SimpleList(lapply(motiflist, to.PWM, type=motif.type, background=background, pseudocount=pseudocount))

  data <- dplyr::bind_rows(by(data, 1:nrow(data), function(x, motifs){
    motif.name <- paste(as.character(x["MotifDB"][[1]]), as.character(x["Motif"][[1]]), as.character(x["provider"][[1]]), sep = "-")
    motif <- motifs[[which(grepl(motif.name, names(motifs)))]]
    x['pvalue.REF'] <- TFMPvalue::TFMsc2pv(motif, as.numeric(x['Ref.score'][[1]]), type = "PWM")
    x['pvalue.ALT'] <- TFMPvalue::TFMsc2pv(motif, as.numeric(x['Alt.score'][[1]]), type = "PWM")
    return(x)
  }, motifs=pwms))
  return(data)
}
