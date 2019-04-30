#' read.input.file
#'
#' Read and strip input file and make ready for motif analysis
#'
#' @param input.file is the path to a file, which can be both a vcf and fasta file.
#' Totally fine if these are zipped.
#' @param ref.genome is the reference genome used to get IRanges and sequence information.
#' @param sample.name which is self explanitory, usefull if you want to codense several
#' outcomes of different samples in a later stage
#' @param ATAC.only can be used to pass a bed file which only has the alterations in an
#' ATAC peak. Based on this bed file, the alterations in input.file will be filtered.
#'
#' @return a GrangesObject with following columns:
#' \item{Sample}{which is the sample.name passed as param}
#' \item{SNP}{name of the SNP}
#' \item{Allel}{Phased allel information if phased vcf is passed}
#' \item{REF}{Reference nucleotide}
#' \item{ALT}{Alternative nucleotide}
#' \item{REF.sequence}{20 nucleotide before and after REF from ref.genome}
#' \item{ALT.sequence}{20 nucleotide before and after ALT from ref.genome}
#'
#' @examples
#' #Load input.files
#' input <- system.file("extdata", "variant.dataset.fasta" , package = "TFBSfindR")
#' #Load reference genome
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' ref.genome <- BSgenome.Hsapiens.UCSC.hg19
#' #Get sample.name from input file
#' sample.name <- "example.dataset"
#' data <- read.input.file(input, ref.genome, sample.name = sample.name, ATAC.only = FALSE)
#'
#' @import Biostrings
#' @import VariantAnnotation
#' @import GenomicRanges
#' @import IRanges
#' @import S4Vectors
#' @import GenomeInfoDb
#' @importFrom SummarizedExperiment rowRanges
#'
#' @export
read.input.file <- function(input.file, ref.genome, sample.name='Unknown', ATAC.only=FALSE){
  if (length(grep(".fasta", input.file)) == 1){
    #library(Biostrings)
    data <- readDNAStringSet(filepath = input.file)
    names <- t(dplyr::bind_cols(lapply(strsplit(names(data), split="|", fixed=TRUE), function(x){return(data.frame(x))})))
    rownames(names) <- 1:length(data)
    names(data) <- names[,6]
    #Should be made so it is not only specific for our fasta files
    data <- GRanges(row.names= 1:length(data),
                    seqnames=names[,3],
                    ranges = IRanges(start=as.numeric(names[,8]),
                                     end=as.numeric(names[,10])),
                    Sample=sample.name,
                    SNP=names[,6],
                    Alteration=names[,12],
                    Ref.Alt=names[,4],
                    sequence=data)
    #add genome data to Granges object
    seqinfo(data) <- attributes(ref.genome)$seqinfo
    data.ref <- data[data$Ref.Alt == "REF"]
    names(mcols(data.ref))[[5]] <- "REF.sequence"
    data.alt <- data[data$Ref.Alt == "ALT"]
    names(mcols(data.alt))[[5]] <- "ALT.sequence"

    #kinda dirty way of not checking snp but just appending the dataframe
    data.ref$ALT.sequence <- data.alt$ALT.sequence
    data <- data.ref[,!(names(mcols(data.ref)) == "Ref.Alt")]
    names(data) <- data$SNP

    if(!(ATAC.only == FALSE)){
      warning("removing snps not in ATAC data")
      OCR <- rtracklayer::import.bed(ATAC.only)
      data <- data[queryHits(findOverlaps(data, OCR))]
      data <- unique(data)
    }
    #To make fasta and vcf the same
    data$Allel <- "*|*"
    Alteration <- strsplit(as.character(data$Alteration), split="")
    data$REF <- do.call(c, lapply(Alteration, function(x) DNAStringSet(x[1])))
    data$ALT <- do.call(c, lapply(Alteration, function(x) DNAStringSet(x[3])))
    data <- data[, c("Sample", "SNP", "Allel", "REF", "ALT", "REF.sequence", "ALT.sequence")]
  }
  if (length(grep(".vcf", input.file)) == 1){
    #library(VariantAnnotation)
    #library(matrixStats)
    genome <- genome(ref.genome)
    data <- readVcf(input.file, genome)
    phased.data <- geno(data)$GT
    data <- rowRanges(data)
    data$Allel <- phased.data
    data$ALT <- unlist(data$ALT)
    data <- keepStandardChromosomes(data)
    data <- data[!(seqnames(data) == "chrX")]
    mcols(data) <- mcols(data)[, c("Allel", "REF", "ALT")]
    data$SNP <- as.character(names(data))
    data$Alteration <- paste(as.character(data$REF), "-", as.character(data$ALT), sep="")
    data$REF.sequence <- getSeq(ref.genome, promoters(data, upstream=20, downstream = 21))
    names(data$REF.sequence) <- names(data)
    #Remove non SNPs
    data <- data[which(nchar(data$Alteration) == 3)]

    if(!(ATAC.only == FALSE)){
      warning("removing snps not in ATAC data")
      OCR <- rtracklayer::import.bed(ATAC.only)
      data <- data[queryHits(findOverlaps(data, OCR))]
      data <- unique(data)
    }

    # NOTE that this is slow -> because of DNAstringSet
    alt <- DNAStringSetList(lapply(data$ALT, function(x) return(x)))
    data$ALT.sequence <- replaceAt(data$REF.sequence, IRanges(21, 21), alt)
    data$Sample <- sample.name
    #To make fasta and vcf the same
    start(ranges(data)) <- start(ranges(data)) - 21
    end(ranges(data)) <- end(ranges(data)) + 20
    data <- data[, c("Sample", "SNP", "Allel", "REF", "ALT", "REF.sequence", "ALT.sequence")]
  }
  #Create minus strand
  strand(data) <- "+"
  data.minus <- data
  data.minus$REF.sequence <- chartr("ACGT", "TGCA", data$REF.sequence)
  data.minus$ALT.sequence <- chartr("ACGT", "TGCA", data$ALT.sequence)
  strand(data.minus) <- "-"
  data <- c(data, data.minus)
  return(data)
}
