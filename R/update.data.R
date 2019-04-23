#' GRanges.update
#'
#' Updates the data to make snp-motif combination readable
#' As well as sorting position at which snp changes motif.
#' Threshold can be used to filter out snp-motif changes
#' that have low scores, making them less relevant.
#'
#' @param data is a single Granges object from \code{TFBS.findR}
#' @param threshold which Delta Kumasaka score should meet.
#' This will make sure that less relavant changes in motif
#' scores will not be reported
#'
#' @return Granges object with columns:
#' \item{Sample}{which is the sample.name passed as param}
#' \item{SNP}{name of the SNP}
#' \item{Allel}{Phased allel information if phased vcf is passed}
#' \item{REF}{Reference nucleotide}
#' \item{ALT}{Alternative nucleotide}
#' \item{Snp.loc}{Location in the motif sequence where variant occurs}
#' \item{Sequence}{Reference sequence matched with motif}
#' \item{MotifDB}{which motif database was used}
#' \item{Motif}{The motif analysed}
#' \item{Ref.score}{PWM compared score for Reference}
#' \item{Alt.score}{PWM compared score for Alternative}
#' \item{Delta.score}{Alt.score - Ref.score}
#' \item{Kuma.ref.score}{posterior probability of transcription
#' factor binding for Reference}
#' \item{Kuma.alt.score}{posterior probability of transcription
#' factor binding for Alternative}
#' \item{Kuma.delta.score}{Kuma.alt.score - Kuma.ref.score}
#'
#' @export
GRanges.update <- function(data, threshold=0.1){
  Delta.Kuma <- data$Kuma.alt.score[[1]] - data$Kuma.ref.score[[1]]
  pwm.length <- 40-length(Delta.Kuma)
  data.c <- c()
  j <- 1
  for(i in 1:length(Delta.Kuma)){
    i.Delta.Kuma <- Delta.Kuma[i]
    if (!(i.Delta.Kuma == 0)){
      i.data <- data
      i.data$Ref.score <- i.data$Ref.score[[1]][i]
      i.data$Alt.score <- i.data$Alt.score[[1]][i]
      i.data$Delta.score <- i.data$Alt.score - i.data$Ref.score
      i.data$Kuma.ref.score <- i.data$Kuma.ref.score[[1]][i]
      i.data$Kuma.alt.score <- i.data$Kuma.alt.score[[1]][i]
      i.data$Kuma.delta.score <- i.Delta.Kuma
      if(i.data$Kuma.delta.score > threshold | i.data$Kuma.delta.score < -threshold){
        i.data$Sequence <- subseq(i.data$REF.sequence, start=i, end=(i+pwm.length))
        start(ranges(i.data)) <- start(ranges(i.data)) + i - 1
        end(ranges(i.data)) <- start(ranges(i.data)) + pwm.length
        i.data$Snp.loc <- 22 - i
        data.c[j] <- list(i.data)
        j <- j + 1
      }
    }
  }
  if(!(is.null(data.c))){
    data <- unlist(GRangesList(data.c))
    data <- data[, c("Sample","SNP","Allel","REF","ALT","Snp.loc","Sequence","MotifDB","Motif","Ref.score","Alt.score","Delta.score","Kuma.ref.score","Kuma.alt.score","Kuma.delta.score")]
    return(data)
  }
}

#' Window.update
#'
#' Update data to a dataframe rather than Granges object
#'
#' @param data is a Granges object from \code{GRanges.update}
#' @param digits number of digits to round off
#'
#' @return dataframe with the following columns:
#' \item{seqnames}{Chromsome}
#' \item{start}{start position of sequence in genome}
#' \item{end}{end position of sequence in genome}
#' \item{width}{size of sequence/motif}
#' \item{strand}{which strand sequence is on}
#' \item{SNP}{variant rs}
#' \item{Allel}{Phased allel information}
#' \item{REF}{Reference nucleotide}
#' \item{ALT}{Alternative nucleotide}
#' \item{Snp.loc}{Location of SNP on sequence}
#' \item{Sequence}{Reference sequence matched with motif}
#' \item{MotifDB}{Motif found in sequence}
#' \item{Ref.score}{PWM compared score for Reference}
#' \item{Alt.score}{PWM compared score for Alternative}
#' \item{Delta.score}{Alt.score - Ref.score}
#' \item{Kuma.ref.score}{posterior probability of transcription
#' factor binding for Reference}
#' \item{Kuma.alt.score}{posterior probability of transcription
#' factor binding for Alternative}
#' \item{Kuma.delta.score}{Kuma.alt.score - Kuma.ref.score}
#'
#' @export
Window.update <- function(data, digits=3){
  data <- as.data.frame(row.names = 1:length(data), data)
  data$seqnames <- as.character(data$seqnames)
  data$strand <- as.character(data$strand)
  numeric.columns <- sapply(data, mode) == 'numeric'
  data[numeric.columns] <- round(data[numeric.columns], digits = digits)
  names(data)[8] <- 'Allel'
  data <- data[order(abs(data$Kuma.delta.score), decreasing=TRUE),]
  rownames(data) <- 1:length(data[,1])
  return(data)
}

#' data.update
#'
#' Updates the data in loop to make snp-motif combination readable
#' As well as sorting position at which snp changes motif.
#' Threshold can be used to filter out snp-motif changes
#' that have low scores, making them less relevant.
#' Lastly, the data will be transformed from a GRanges object
#' to a dataframe.
#'
#' @param data is a single Granges object from \code{TFBS.findR}
#' @param threshold which Delta Kumasaka score should meet.
#' This will make sure that less relavant changes in motif
#' scores will not be reported
#'
#' @return dataframe with the following columns:
#' \item{seqnames}{Chromsome}
#' \item{start}{start position of sequence in genome}
#' \item{end}{end position of sequence in genome}
#' \item{width}{size of sequence/motif}
#' \item{strand}{which strand sequence is on}
#' \item{SNP}{variant rs}
#' \item{Allel}{Phased allel information}
#' \item{REF}{Reference nucleotide}
#' \item{ALT}{Alternative nucleotide}
#' \item{Snp.loc}{Location of SNP on sequence}
#' \item{Sequence}{Reference sequence matched with motif}
#' \item{MotifDB}{Motif found in sequence}
#' \item{Ref.score}{PWM compared score for Reference}
#' \item{Alt.score}{PWM compared score for Alternative}
#' \item{Delta.score}{Alt.score - Ref.score}
#' \item{Kuma.ref.score}{posterior probability of transcription
#' factor binding for Reference}
#' \item{Kuma.alt.score}{posterior probability of transcription
#' factor binding for Alternative}
#' \item{Kuma.delta.score}{Kuma.alt.score - Kuma.ref.score}
#'
#' @export
data.update <- function(data, threshold=0.1){
  #Update the GrangesList object
  data <- unlist(GRangesList(unlist(lapply(data, GRanges.update, threshold=threshold), use.names = FALSE)), use.names = FALSE)
  data <- data[order(data$Kuma.delta.score, decreasing=TRUE)]
  data <- Window.update(data)
  return(data)
}
