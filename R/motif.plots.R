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
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @export
snp.plot <- function(snp, method='both', motif="MEF2A", strand=c("-", "+")){
  #library(ggplot2)
  #library(cowplot)

  #Because CMD check
  Alt.score <- Ref.score <- Delta <- Kuma.Alt.score <- Kuma.Ref.score <- Kuma.Delta <- NULL

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
  if (method=="Kuma" | method =="both"){
    snp.frame$Kuma.Ref.score <- snp$Kuma.ref.score[[1]]
    snp.frame$Kuma.Alt.score <- snp$Kuma.alt.score[[1]]
    snp.frame$Kuma.Delta <- snp$Kuma.alt.score[[1]]-snp$Kuma.ref.score[[1]]
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
  #+ geom_line(mapping = aes(y=Kuma.Alt.score, colour="red")) + geom_segment(aes(x=max.at.Kuma, y= -10, xend=max.at.Kuma, yend=10), colour='black', linetype="dashed", size=0.1)
  plot.1 <- ggplot(snp, mapping = aes(x=1:nrow(snp))) + geom_line(mapping = aes(y=Alt.score, colour="Alt")) + geom_line(mapping = aes(y=Ref.score, colour="Ref")) + labs(x="motif index", y="PWM score")
  plot.2 <- ggplot(snp, mapping = aes(x=1:nrow(snp), y=Delta)) + geom_line() + labs(x="motif index", y=expression(Delta ~ "PWM score"))
  if (toupper(method)=="KUMA" | toupper(method) =="BOTH"){
    if (max(snp$Kuma.Ref.score) > max(snp$Kuma.Alt.score)){
      max.at.Kuma <- which(snp$Kuma.Ref.score == max(snp$Kuma.Ref.score))
    } else{
      max.at.Kuma <- which(snp$Kuma.Alt.score == max(snp$Kuma.Alt.score))
    }
    plot.3 <- ggplot(snp, mapping = aes(x=1:nrow(snp))) + geom_line(mapping = aes(y=Kuma.Alt.score, colour="Alt"))  + geom_line(mapping = aes(y=Kuma.Ref.score, colour="Ref")) + labs(x="motif index", y="Kuma score")
    plot.4 <- ggplot(snp, mapping = aes(x=1:nrow(snp), y=Kuma.Delta)) + geom_line() + labs(x="motif index", y=expression(Delta ~ "Kuma score"))
    return(plot_grid(plot.1, plot.2, plot.3, plot.4, labels = "auto", align = "h"))
  }
  return(plot_grid(plot.1, plot.2, labels = "auto", align = "h"))
}

#' variant.plot
#'
#' Plots for a specific snp, the PWM scores and Kumasaka of all motifs in the variant
#'
#' @param snp is a specific item from TFBS.findR after data.update which include only
#' one snp.
#' @param MotifDb is a simplelist of all the motifs PPM provided within snp.
#' @param method is the method used for TFBs.findR. if both or Kuma is provided
#' kumasaka scores will also be plotted.
#'
#' @return one or two plot for PWM score and/or Kuma score for a specific snp.
#'
#' @import ggplot2
#' @import GenomicRanges
#' @import gridExtra
#' @import ggrepel
#' @export
variant.plot <- function(snp, MotifDb, method="both"){
  #Change to Granges object
  if(toupper(method)=="BOTH" | toupper(method)=="KUMA"){
    snp <- GRanges(seqnames = snp$seqnames, ranges=IRanges(start=snp$start, end=snp$end, width = snp$width), strand=snp$strand, SNP=snp$SNP, REF=snp$REF, ALT=snp$ALT, Snp.loc=snp$Snp.loc, Sequence=snp$Sequence, MotifDB=snp$MotifDB, provider=snp$provider, Motif=snp$Motif, Ref.score=snp$Ref.score, Alt.score=snp$Alt.score, Delta.score=snp$Delta.score, Kuma.ref.score=snp$Kuma.ref.score, Kuma.delta.score=snp$Delta.score)
  } else{
    snp <- GRanges(seqnames = snp$seqnames, ranges=IRanges(start=snp$start, end=snp$end, width = snp$width), strand=snp$strand, SNP=snp$SNP, REF=snp$REF, ALT=snp$ALT, Snp.loc=snp$Snp.loc, Sequence=snp$Sequence, MotifDB=snp$MotifDB, provider=snp$provider, Motif=snp$Motif, Ref.score=snp$Ref.score, Alt.score=snp$Alt.score, Delta.score=snp$Delta.score)
  }

  #Get high impact
  snp <- unlist(GRangesList(lapply(snp, function(x, y){
    motif <- y[which(grepl(as.character(x$Motif), names(y)))]
    if(length(motif) > 1){
      index <- as.numeric(which(unlist(strsplit(names(motif), split="-")) == as.character(x$Motif)) < 5)
      if(index == 0){index<-2}
      motif <- motif[index]
    }
    motif <- motif[[1]][,x$Snp.loc]
    ref <- motif[[as.numeric(chartr("ACGT", "1234", as.character(x$REF)))]]
    alt <- motif[[as.numeric(chartr("ACGT", "1234", as.character(x$ALT)))]]
    if(abs(alt - ref) > 0.75){
      x$high.impact <- 'yes'
    } else{
      x$high.impact <- 'no'
    }
    return(x)
  }, y=MotifDb)), use.names = FALSE)

  snp <- as.data.frame(snp)

  p1 <- ggplot(snp, mapping=aes(xmin=start-0.5, xmax=end+0.5, ymin=Delta.score-0.1, ymax=Delta.score+0.1, label=Motif)) + geom_vline(xintercept = snp$start+snp$Snp.loc, alpha=0.1) + geom_rect(alpha=0, colour="black") + scale_x_continuous(labels = c("Start", "SNP", "End"), breaks = c((min(snp$start)+max(snp$start))/2-25, (min(snp$start)+max(snp$start))/2, (min(snp$start)+max(snp$start))/2+25), limits = c((min(snp$start)+max(snp$start))/2-25, (min(snp$start)+max(snp$start))/2+25)) + geom_text_repel(data = snp[snp$Snp.loc > 7,], aes(x=end+0.5, y=Delta.score), nudge_x = 20, segment.size  = 0.2, segment.color = "grey50", direction = "y", hjust = 0) + geom_text_repel(data = snp[snp$Snp.loc <= 7,], aes(x=start-0.5, y=Delta.score), nudge_x = -20, segment.size = 0.2, segment.color = "grey50", direction = "y", hjust = 0) + theme(axis.ticks = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.1), axis.title.x = element_blank()) + ylab(expression(Delta ~ "position weight matrix change"))
  if(any(snp$high.impact == "yes")){
    p1 <- p1 + geom_rect(data=snp[snp$high.impact == "yes",], aes(xmin=start+Snp.loc-0.5, xmax=start+Snp.loc+0.5, ymin=Delta.score-0.1, ymax=Delta.score+0.1), colour="black")
  }

  if(toupper(method)=="BOTH" | toupper(method)=="KUMA"){
    p2 <- ggplot(snp, mapping=aes(xmin=start-0.5, xmax=end+0.5, ymin=Kuma.delta.score-0.1, ymax=Kuma.delta.score+0.1, label=Motif)) + geom_vline(xintercept = snp$start+snp$Snp.loc, alpha=0.1) + geom_rect(alpha=0, colour="black") + scale_x_continuous(labels = c("Start", "SNP", "End"), breaks = c((min(snp$start)+max(snp$start))/2-25, (min(snp$start)+max(snp$start))/2, (min(snp$start)+max(snp$start))/2+25), limits = c((min(snp$start)+max(snp$start))/2-25, (min(snp$start)+max(snp$start))/2+25)) + geom_text_repel(data = snp[snp$Snp.loc > 7,], aes(x=end, y=Kuma.delta.score), nudge_x = 20, segment.size  = 0.2, segment.color = "grey50", direction = "y", hjust = 0) + geom_text_repel(data = snp[snp$Snp.loc <= 7,], aes(x=start, y=Kuma.delta.score), nudge_x = -20, segment.size = 0.2, segment.color = "grey50", direction = "y", hjust = 0) + theme(axis.ticks = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.1), axis.title.x = element_blank()) + ylab(expression(Delta ~ "Kumasaka score change"))
    if(any(snp$high.impact == "yes")){
      p2 <- p2 + geom_rect(data=snp[snp$high.impact == "yes",], aes(xmin=start+Snp.loc-0.5, xmax=start+Snp.loc+0.5, ymin=Kuma.delta.score-0.1, ymax=Kuma.delta.score+0.1), colour="black")
    }
    return(grid.arrange(p1, p2, ncol=2))
  }
  return(p1)
}
