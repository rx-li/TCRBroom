#' Title Plot heatmap for pairwise TCR overlaps
#'
#' @description The function is used to plot heatmap for pairwise TCR overlaps among
#' the subjects. Each subject can have multiple samples. It may take a while if
#' data is big.
#'
#' @param meta  A meta table that contains all samples' names and the corresponding
#' patients' ids. The first column should be sample name, and second column should
#' be subjects' id
#' @param dir Data path
#'
#' @import RColorBrewer gplots
#'
#' @return A heatmap and a table presenting pairwise TCR overlaps between any 2 subjects.
#'
#' @export
#'
#' @examples # pairwise_overlap <- checkSeq(meta, "./")
#'
source("utils.R")
checkPairwise <- function(meta, dir){
  uniqID <- unique(meta[, 2])
  overlapmat <- matrix(NA, length(uniqID), length(uniqID))
  for(jj in 1:(length(uniqID)-1)) {
    for(qq in (jj+1):length(uniqID)) {
      onesum1 <- GetOneSumTable(uniqID[jj], meta=meta, path=dir)
      onesum2 <- GetOneSumTable(uniqID[qq], meta=meta, path=dir)

      overlapmat[jj, qq] <- length(intersect(rownames(onesum1), rownames(onesum2)))
      overlapmat[jj, jj] <- nrow(onesum1)
      overlapmat[qq, qq] <- nrow(onesum2)
    }
  }
  rownames(overlapmat) = colnames(overlapmat) = uniqID
  overlapmat2 <- overlapmat
  for(i in seq(nrow(overlapmat))){
    for (j in seq(nrow(overlapmat))){
      if (i < j){
        overlapmat2[i, j] <- overlapmat[i, j] / (overlapmat[i, i] +
                                                   overlapmat[j, j] - overlapmat[i, j])
      }
    }
  }
  diag(overlapmat2) <- 1
  propmat <- signif(overlapmat2, 2)
  diag(propmat) <- 0
  Colors=rev(brewer.pal(11,"Spectral"))
  Colors=colorRampPalette(Colors)(100)
  heatmap.2(propmat, Rowv=FALSE, Colv=FALSE, dendrogram='none', trace="none",
            col=Colors)
  return(propmat)
}
