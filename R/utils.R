#' Title Read data in tsv format for one patient.
#'
#' @description Each patient can have multiple samples. The data should have 3 columns:
#' sequences, sequence counts and sequences frequencies.
#'
#' @param thismeta A meta table contains the name of the samples for the patient.
#' @param data_dir Folder location for data.
#'
#' @return A dataframe that contains sequences, sequence counts and sequence frequencies.
#'
ReadOnePat <- function(thismeta, data_dir) {
  onealldata <- list()
  for(i in 1:nrow(thismeta)) {
    onealldata[[i]] <- read.table(paste0(data_dir, thismeta[i,1]),
                                  colClasses = c("character", "numeric", "numeric"),
                                  header = TRUE, sep = "\t")
  }
  return(list(onealldata = onealldata,
              thismeta = thismeta))
}


#' Title Summary all sequences in a table for a patient
#'
#' @description Each patient can have multiple samples. The data should have 3 columns:
#' sequences, sequence counts and sequences frequencies.
#'
#' @param thisID Patient ID
#' @param meta A meta table that contains all patients' id and the corresponding
#' samples' names. The first column should be sample names, and second column should
#' be subjects' ids.
#'
#' @return A dataframe that contains sequences, sequence counts and sequence frequencies.
#'
GetOneSumTable <- function(thisID, meta, path) {
  thismeta <- meta[meta[, 2] == thisID, ]
  oneout <- ReadOnePat(thismeta, path)

  allseq <- c()
  for(i in 1:nrow(thismeta)) {
    allseq <- c(allseq, oneout$onealldata[[i]][,1])
  }
  #length(allseq)
  uniqLT <- unique(allseq)
  counttable <- matrix(0, length(uniqLT), nrow(thismeta)*2)
  for(i in 1:nrow(thismeta)) {
    tmp <- oneout$onealldata[[i]]
    idx0 <- match(uniqLT, tmp$nucleotide)
    counttable[which(!is.na(idx0)), c(2*(i-1)+1, 2*(i-1)+2)] <- as.matrix(tmp[na.omit(idx0), 2:3])
  }
  rownames(counttable) <- uniqLT
  colnames(counttable) <- paste0(rep(oneout$thismeta$sample_id, each = 2), c("_count", "_freq"))
  return(counttable)
}
