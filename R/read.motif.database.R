#' read.motif.database.R
#'
#' Read and strip text motif database files for motif analysis
#'
#' @param file is a path to a .txt file with motif scores
#'
#' @return a SimpleList with motif matrixes
#'
#' @export
read.motif.database <- function(file){
  data <- readLines(file)
  index <- grep(">", data)
  motifs <- c()
  for (i in 1:length(index)){
    name <- data[index[i]]
    #PROBLEM WITH LAST 1
    if (is.na(unlist(index[i+1]))){
      second.index <- length(data)
    } else{
      second.index <- index[i+1]-1
    }
    data.i <- data[(index[i]+1):second.index]
    length.i <- length(data.i)
    data.i <- as.numeric(unlist(strsplit(data.i, split="\t")))
    data.i <- matrix(data=data.i, ncol=length.i, nrow=4)
    colnames(data.i) <- 1:length.i
    rownames(data.i) <- list("A", "C", "G", "T")
    motifs[name] <- list(data.i)
  }
  data <- SimpleList(motifs)
  return(data)
}
