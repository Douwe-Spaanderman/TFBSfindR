#' snp.plot
#'
#' Plots for a specific snp, strand and motif the PWM scores and Kumasaka
#' scores if provided for Ref vs Alt
#'
#' @param snp is a specific GRanges item from TFBS.findR which include only
#' one snp.
#' @param method is the method used for TFBs.findR. if both or Kuma is provided
#' kumasaka scores will also be plotted.
#' @param motif is the motif which should be investigated from the Granges object.
#' Note that the motif needs to be present in the GRanges object.
#' @param strand is the strand which should be investigated from the Granges
#' Object. As default, both will be accepted, however if both plus and minus
#' strand are present, than ploting doens't work.
#'
#' @return two or four plots with PWM scores, delta PWM scores, Kumasaka
#' scores and delta Kumasaka scores
#'
#' @examples
#' #import data from TFBS.findR
#' snp <- data[data$SNP %in% "rs74813460"]
#' snp.plot(snp, method="both", motif="MYC::MAX", strand="+")
#'
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @export
snp.plot <- function(snp, method='both', motif="MEF2A", strand=c("-", "+")){
  library(ggplot2)
  library(cowplot)

  if(length(motif) > 1){
    stop("please provide only one motif to analyse")
  }
  snp <- snp[snp$Motif %in% motif]
  snp <- snp[strand(snp) %in% strand]

  if(length(snp) > 1){
    stop("please also provide a specific strand as both strand have snp-motif combination")
  }

  snp.frame <- data.frame(name=snp$Motif,
                          Ref.score=snp$Ref.score[[1]],
                          Alt.score=snp$Alt.score[[1]],
                          Delta=snp$Alt.score[[1]]-snp$Ref.score[[1]])
  if (method=="gaf" | method =="both"){
    snp.frame$Gaf.Ref.score <- snp$Gaf.ref.score[[1]]
    snp.frame$Gaf.Alt.score <- snp$Gaf.alt.score[[1]]
    snp.frame$Gaf.Delta <- snp$Gaf.alt.score[[1]]-snp$Gaf.ref.score[[1]]
  }
  snp <- snp.frame

  if(all(!(is.na(snp))) == FALSE){
    warning("NaN values are transformed to 0")
    snp[is.na(snp)] <- 0
  }

  if (max(snp$Ref.score) > max(snp$Alt.score)){
    max.at <- which(snp$Ref.score == max(snp$Ref.score))
  } else{
    max.at <- which(snp$Alt.score == max(snp$Alt.score))
  }
  #+ geom_line(mapping = aes(y=Gaf.Alt.score, colour="red")) + geom_segment(aes(x=max.at.gaf, y= -10, xend=max.at.gaf, yend=10), colour='black', linetype="dashed", size=0.1)
  plot.1 <- ggplot(snp, mapping = aes(x=1:nrow(snp))) + geom_line(mapping = aes(y=Alt.score, colour="Alt")) + geom_line(mapping = aes(y=Ref.score, colour="Ref")) + labs(x="motif index", y="PWM.score")
  plot.2 <- ggplot(snp, mapping = aes(x=1:nrow(snp), y=Delta)) + geom_line() + labs(x="motif index", y="Delta.PWM.score")
  if (toupper(method)=="KUMA" | toupper(method) =="BOTH"){
    if (max(snp$Gaf.Ref.score) > max(snp$Gaf.Alt.score)){
      max.at.gaf <- which(snp$Gaf.Ref.score == max(snp$Gaf.Ref.score))
    } else{
      max.at.gaf <- which(snp$Gaf.Alt.score == max(snp$Gaf.Alt.score))
    }
    plot.3 <- ggplot(snp, mapping = aes(x=1:nrow(snp))) + geom_line(mapping = aes(y=Gaf.Alt.score, colour="Alt"))  + geom_line(mapping = aes(y=Gaf.Ref.score, colour="Ref")) + labs(x="motif index", y="Gaf.score")
    plot.4 <- ggplot(snp, mapping = aes(x=1:nrow(snp), y=Gaf.Delta)) + geom_line() + labs(x="motif index", y="Delta.Gaf.score")
    return(plot_grid(plot.1, plot.2, plot.3, plot.4, labels = "auto", align = "h"))
  }
  return(plot_grid(plot.1, plot.2, labels = "auto", align = "h"))
}
