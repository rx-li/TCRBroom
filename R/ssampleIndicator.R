#' Title Small sample indicator  
#'
#' @description The function is used to identify if a subject has small number of 
#' frequencies. 
#'
#' @param meta  A meta table that contains all samples' names and the corresponding
#' patients' ids. The first column should be sample name, and second column should
#' be subjects' id
#' @param dir Data path
#' @param low.seq A integer that define "low frequency" for each subject. 
#' If a subject has sequences lower than or equal to this number, it will be noted as "Ouch!"; 
#' Otherwise it will show NA. 
#'
#' @import RColorBrewer gplots
#'
#'
#' @export
#'
#' @examples # 
#'
ssampleIndicator <- function(meta, dir, low.seq=1000){
  uniqID <- unique(meta[, 2])
  uniqSeq <- list()
  for(j in seq(length(uniqID))) {
    onesum <- GetOneSumTable(uniqID[j], meta=meta, path=dir)
    uniqSeq[[j]] <- nrow(onesum)
  }
  d <- data.frame(subject_id=uniqID, 
                  n_unique_seq=unlist(uniqSeq))
  d$low_seq <- ifelse(d$n_unique_seq <= low.seq, "Alert", NA)
  return(d)
}
