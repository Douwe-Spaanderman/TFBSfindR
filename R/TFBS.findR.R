#' to.PWM
#'
#' Function to change position frequency matrix/position count matrix or
#' position probability matrix to a position weight matrix.
#'
#' Several methods have been discussed to calculate PWM,
#' Papatsenko et al., 2003, suggested the method we use now, with the use of
#' a pseudocount and sum log. As suggested by papatsenko, pseudocount is
#' calculated by log(sequenceCount), in which sequenceCount is the amount
#' of times PFM is sequenced. If however this information isn't available
#' -or the other method is prefered- than the pseudocount can be set
#' beforehand. Nakai et al., 2008, discussed than smaller pseudocounts
#' preserve data model better than larger pseudocount. The lowest value
#' they analysed was 0.01 and therefore this is also our default.
#' Note that if you want to calculate the pseudocount based on the log of
#' sequence reads, this param should be set to log.of.reads.
#'
#' @param x is a GrangesObject from \code{read.input.file}
#' @param type is the type of input file, which can be "PFM", "PCM",
#' "PPM" and "PWM". If the latter is passed, than the PWM is returned
#' without changes
#' @param background is the background probability for each nucleotide.
#' @param pseudocount value to change the outcome model, which removes
#' possibility for -inf values.
#'
#' @return transformed matrix of x based on the method selected
#'
#' @export
to.PWM <- function(x, type="PFM", background=c(A=0.25, C=0.25, G=0.25, T=0.25), pseudocount=0.01){
  #PWM calculations are done as https://academic.oup.com/nar/article/31/20/6016/1039515
  #Pseudocount is set to 0.01 because of https://academic.oup.com/nar/article/37/3/939/1074640
  #Pseudocount can also be calculated with log(sequenceCount) using log.of.reads for pseudocount
  pwm.name <- names(x)[1]

  if(toupper(type) == "PFM" | toupper(type) == "PCM"){
    if(pseudocount == "log.of.reads"){
      sequenceCount <- sum(x[,1])
      pseudocount <- log(sequenceCount)
      pwm <- log((x+pseudocount*background)/((sequenceCount+pseudocount)*background))
    } else{
      x <- x/sum(x[,1])
      pwm <- log((x + (pseudocount*background))/((1 + pseudocount)*background))
    }
  }
  if(toupper(type) == "PPM"){
    if(pseudocount == "log.of.reads"){
      sequenceCount <- as.numeric(names(x)[2])
      x <- x*sequenceCount
      pseudocount <- log(sequenceCount)
      pwm <- log((x+pseudocount*background)/((sequenceCount+pseudocount)*background))
    } else{
      pwm <- log((x+(pseudocount*background))/((1+pseudocount)*background))
    }
  }
  if(toupper(type) == "PWM"){
    warning('no changes have been made to motiflist as "PWM" was provided')
    pwm <- x
  }

  names(pwm) <- pwm.name
  return(pwm)
}

#' pwm.compare
#'
#' Compare position weight matrix with both the reference
#' and alternative sequence in data, in order to calculate
#' score for each position
#'
#' Here we analyse our reference and alternative sequence
#' based on a 'rolling TF' model. Calculated pwm from
#' \code{to.PWM} is compared to Grangesobject from
#' \code{read.input.file}. PWM are compared to every possible
#' position in the sequence, even outside SNP.
#'
#' @param pwm are calculated pwm scores from \code{to.PWM}
#' @param data are GRanges objects from \code{read.input.file}
#'
#' @return GrangesObject with added columns:
#' \item{Ref.score}{list of PWM compared score for each position
#' to Ref.sequence}
#' \item{Alt.score}{list of PWM compared score for each position
#' to Alt.sequence}
#' Note that the amount of GrangesObject are also *amount of pwms,
#' as a score list is provided for each individual pwm
#'
#' @export
pwm.compare <- function(pwm, data){
  pwm.name <- names(pwm)[1]
  pwm.name <- strsplit(pwm.name, split="-", fixed=TRUE)[[1]]
  names(pwm) <- NULL

  #ref seq
  ref <- data$REF.sequence
  ref <- strsplit(as.character(ref), split="")[[1]]
  ref <- as.numeric(chartr("ACGT", "1234", ref))

  #alt seq
  alt <- data$ALT.sequence
  alt <- strsplit(as.character(alt), split="")[[1]]
  alt <- as.numeric(chartr("ACGT", "1234", alt))

  if (as.character(strand(data)) == "-"){
    ref <- rev(ref)
    alt <- rev(alt)
  }

  length.motif <- length(colnames(pwm))

  pwm.ref <- c()
  pwm.alt <- c()
  for (i in 1:(length(ref)-(length.motif))){
    ref.i <- ref[i:(i+length.motif)]
    alt.i <- alt[i:(i+length.motif)]

    pwm.ref.i <- c()
    pwm.alt.i <- c()
    for(j in 1:length.motif){
      pwm.ref.i[j] <- pwm[ref.i[j],j]
      pwm.alt.i[j] <- pwm[alt.i[j],j]
    }

    pwm.ref[i] <- sum(pwm.ref.i)
    pwm.alt[i] <- sum(pwm.alt.i)
  }

  data$MotifDB <- pwm.name[2]
  data$provider <- pwm.name[4]
  data$Motif <- pwm.name[3]
  data$Ref.score <- list(pwm.ref)
  data$Alt.score <- list(pwm.alt)
  return(data)
}

#' analyse.pwm
#'
#' Analyse position weight matrix to Grangesobject from \code{read.input.file}.
#' This is only to do analyses in loop for each motif in motiflist.
#'
#' @param data is a single Grangesobject from \code{read.input.file}
#' @param motiflist is a simplelist which include several (at least 2)
#' motifs which are position weight matrices, position probabilty matrices or
#' position frequency matrices.
#' @param motif.type is the type of matrices passed in motiflist. Can be "PFM",
#' "PCM", "PPM", "PWM".
#' @param background is the background probability for each nucleotide.
#' @param pseudocount value to change the outcome model, which removes
#' possibility for -inf values.
#'
#' @return GrangesObject with added columns:
#' \item{Ref.score}{list of PWM compared score for each position
#' to Ref.sequence}
#' \item{Alt.score}{list of PWM compared score for each position
#' to Alt.sequence}
#' Note that the amount of GrangesObject are also *amount of pwms,
#' as a score list is provided for each individual pwm
#'
#' @import S4Vectors
#' @export
analyse.pwm <- function(data, motiflist, motif.type="PFM", background=c(A=0.25, C=0.25, G=0.25, T=0.25), pseudocount=0.01){
  #Add names to pwms
  for (i in 1:length(motiflist)){
    names(motiflist[[i]]) <- names(motiflist[i])
    names(motiflist[[i]])[2] <- mcols(motiflist)$sequenceCount[i]
  }
  #Changing motifs to PWMs
  pwms <- SimpleList(lapply(motiflist, to.PWM, type=motif.type, background=background, pseudocount=pseudocount))

  #Comparing sequence to PMWs
  data <- unlist(GRangesList(lapply(pwms, pwm.compare, data=data)), use.names = FALSE)
}

#' kumasaka.score
#'
#' kumaska et al., 2019, used a different method of scoring
#' in which posterior probability of transcription factor
#' binding is calculated
#' MORE INFO
#'
#' @param data is a Grangesobject from \code{analyse.pwm}
#' @param background is the postion frequency matrix with each
#' background score of a nucleotide
#' @param prior is set the False discovery rate (FDR) which is
#' set to 10%. 0.1 is the default setting as we normally pass
#' a dataset from 1 sample.
#'
#' @return Grangesobject with added columns:
#' \item{Kuma.ref.score}{list of posterior probability of transcription
#' factor binding score for each position and motif to Ref.sequence}
#' \item{Kuma.alt.score}{list of posterior probability of transcription
#' factor binding score for each position and motif to Alt.sequence}
#'
#' @export
kumasaka.score <- function(data, background=c(A=0.25, C=0.25, G=0.25, T=0.25), prior=0.1){
  pwm.length <- nchar(data$REF.sequence) - length(data$Ref.score[[1]])

  #Check if names consists of only ACGT and get them in that order -> this is to correct for different inputs
  if(length(names(background)) > 0){
    if(all(names(background) %in% c("A", "C", "G", "T")) == TRUE & Reduce("+",nchar(names(background))) == 4){
      background <- background[c("A", "C", "G", "T")]
    } else{
      background = c(A=0.25, C=0.25, G=0.25, T=0.25)
    }
  } else{
    names(background) <- c("A", "C", "G", "T")
  }

  if (all(background == c(A=0.25, C=0.25, G=0.25, T=0.25))){
    background.score <- 0.25*pwm.length
  } else{
    background.PWM <- matrix(data = background, nrow = 4, ncol = pwm.length, dimnames = list(c("A", "C", "G", "T"), c(1:pwm.length)))
    ### CURRENTLY NOT WORKING BACKGROUND SCORE IF OTHER BACKGROUND THAN 0.25
    warning('currently other background not implemented')
  }

  #Have to loop through because negative values should return 0
  Ref.score <- c()
  Alt.score <- c()
  for(i in 1:length(data$Ref.score[[1]])){
    x.Ref <- data$Ref.score[[1]][i]
    x.Alt <- data$Alt.score[[1]][i]

    if (x.Ref < 0){
      Ref.score[i] <- 0
    } else{
      Ref.score[i] <- (prior*x.Ref)/(((1-prior)*background.score)+(prior*x.Ref))
    }

    if (x.Alt < 0){
      Alt.score[i] <- 0
    } else{
      Alt.score[i] <- (prior*x.Alt)/(((1-prior)*background.score)+(prior*x.Alt))
    }
  }

  data$Kuma.ref.score <- list(Ref.score)
  data$Kuma.alt.score <- list(Alt.score)
  return(data)
}

#' TFBS.findR
#'
#' TFBS.findR analyses Grangesobjects from \code{read.input.file}
#' if they match motifs in the motiflist.
#'
#' Motif-sequence analyses is conducted using several steps:
#' \code{analyse.pwm}. First changes motifs provided in motiflist
#' to position weight matrices. Secondly, these PWMs are compared
#' to individual Grangesobjects from \code{read.input.file} and
#' PWM scores are given for each position in the sequence of these
#' Grangesobjects
#' \code{kumasaka.score} analyses values provided by analyse.pwm
#' and calculates the posterior probability of transcription factor
#' binding. These values are than similarly passed as a list to the
#' Granges object
#'
#' @param data is a GRangesobject from \code{read.input.file}
#' @param motiflist is a simplelist of motifs
#' @param motif.type is the type of motifs passed. These can be
#' "PWM", "PFM", "PCM" or "PPM".
#' @param method can be used to determine if kumasaka score will be
#' calculated.
#' @param background is the background probability for each nucleotide.
#' @param pseudocount value to change the outcome model, which removes
#' possibility for -inf values.
#' @param prior is set the False discovery rate (FDR) which is
#' set to 10%. 0.1 is the default setting as we normally pass
#' a dataset from 1 sample.
#' @param BPPARAM are settings by BiocParallel, which can be used
#' to do multiple core calculations. bpparam() are the current settings
#' on your system.
#'
#' @return a Grangesobject with following columns:
#' \item{Sample}{which is the sample.name passed as param}
#' \item{SNP}{name of the SNP}
#' \item{Allel}{Phased allel information if phased vcf is passed}
#' \item{REF}{Reference nucleotide}
#' \item{ALT}{Alternative nucleotide}
#' \item{REF.sequence}{20 nucleotide before and after REF from ref.genome}
#' \item{ALT.sequence}{20 nucleotide before and after ALT from ref.genome}
#' \item{MotifDB}{which motif database was used}
#' \item{provider}{provider from the motif}
#' \item{Motif}{The motif analysed}
#' \item{Ref.score}{list of PWM compared score for each position
#' to Ref.sequence}
#' \item{Alt.score}{list of PWM compared score for each position
#' to Alt.sequence}
#' \item{Kuma.ref.score}{list of posterior probability of transcription
#' factor binding score for each position and motif to Ref.sequence}
#' \item{Kuma.alt.score}{list of posterior probability of transcription
#' factor binding score for each position and motif to Alt.sequence}
#'
#' @import BiocParallel
#' @export
TFBS.findR <- function(data, motiflist, motif.type="PFM", method="both", background=c(A=0.25, C=0.25, G=0.25, T=0.25), pseudocount=0.01, prior=0.1, BPPARAM=bpparam()){
  #Analyse pwms on occurence of motifs
  data <- unlist(GRangesList(bplapply(data, analyse.pwm, motiflist=motiflist,
                                      motif.type=motif.type, background=background,
                                      pseudocount=pseudocount, BPPARAM=BPPARAM)), use.names = FALSE)
  if (toupper(method)=="KUMA" | toupper(method) =="BOTH"){
    #Calculating Kumasaka score
    data <- unlist(GRangesList(lapply(data, kumasaka.score, background=background, prior=prior)), use.names=FALSE)
  }
  return(data)
}
