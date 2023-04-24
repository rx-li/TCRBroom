#' Check overlaps between two samples
#' 
#' Given two subject's names, the function will check the overlaps between all their
#' samples. 
#'
#' @param meta Meta table 
#' @param id1 The first subject
#' @param id2 The second subject
#'
#' @return A matrix that contains 
#' @export
#'
#' @examples
#' # check2Samples(meta, id1="P01", id2="P02")
#' 
check2Samples <- function(meta, id1, id2, data_dir) {
  
  thismeta <- meta[meta$subject_id == id1, ]
  oneout1 <- ReadOnePat(thismeta, data_dir)
  
  thismeta <- meta[meta$subject_id == id2, ]
  oneout2 <- ReadOnePat(thismeta, data_dir)
  
  mysummary <- matrix(NA, length(oneout1[[1]])+1, length(oneout2[[1]])+1)
  for(i in 1:length(oneout1[[1]])) {
    for(j in 1:length(oneout2[[1]])) {
      mysummary[i,j] <- length(intersect(oneout1[[1]][[i]]$nucleotide, oneout2[[1]][[j]]$nucleotide))
      mysummary[i,j+1] <- length(oneout1[[1]][[i]]$nucleotide)
      mysummary[i+1,j] <- length(oneout2[[1]][[j]]$nucleotide)
    }
  }
  rownames(mysummary) <- c(oneout1$thismeta$sample_name, "sampleTotal")
  colnames(mysummary) <- c(oneout2$thismeta$sample_name, "sampleTotal")
  
  # calculate the pct
  mysummary_pct <- as.matrix(mysummary)
  for(i in seq(nrow(mysummary_pct))){
    for (j in seq(ncol(mysummary_pct))){
      mysummary_pct[i, j] <- round(
        mysummary_pct[i, j] / (mysummary_pct[i, ncol(mysummary_pct)] + 
                                 mysummary_pct[nrow(mysummary_pct), j] - 
                                 mysummary_pct[i, j]), 
        5)
    }
  }
  mysummary_pct[nrow(mysummary_pct), ] <- mysummary[nrow(mysummary), ]
  mysummary_pct[, ncol(mysummary_pct)] <- mysummary[, ncol(mysummary)]
  
  return(list(overlap_n=mysummary, overlap_pct=mysummary_pct))
}