#' Title Generate oncoprint for cross-cohort overlap
#'
#' @param meta A meta table that contains all samples' names and the corresponding
#' patients' ids. The first column should be sample name, and second column should
#' be subjects' id
#' @param dir Data path
#'
#' @import ComplexHeatmap
#'
#' @return A OncoPrint and a table presenting sequence prevalence rate among the subjects.
#'
#' @export
#'
#' @examples # crossc_overlap <- checkSeq(meta, "./")
checkCrosscohort <- function(meta, dir){
  recordSeq <- c()
  for(i in 1:nrow(meta)) {
    oneseq <- read.table(paste0(dir, meta[i,1]),
                         colClasses = c("character", "numeric", "numeric"),
                         header = TRUE, sep = "\t")
    recordSeq <- unique(c(recordSeq, oneseq$nucleotide))
  }

  sumTable <- matrix(NA, length(recordSeq), length(unique(meta$subject_id)))
  uniqID <- unique(meta$subject_id)
  for(i in 1:length(uniqID)) {

    onetab <- GetOneSumTable(thisID = unique(meta[, 2])[i], meta, dir)

    idx11 <- match(recordSeq, rownames(onetab))
    sumTable[!is.na(idx11), i] <- 1
  }
  rownames(sumTable) <- recordSeq
  colnames(sumTable) <- unique(meta$subject_id)
  prev_sum <- rowSums(sumTable,na.rm = TRUE)
  prev_prop <- rowSums(sumTable,na.rm = TRUE)/length(uniqID)
  subTable3 <- cbind(sumTable, prev_sum, prev_prop)
  sumTable <- sumTable[order(subTable3[, length(uniqID)+1], decreasing = TRUE), ]
  subTable3 <- subTable3[order(subTable3[, length(uniqID)+1], decreasing = TRUE), ]

  sumTable[is.na(sumTable)] = ""
  sumTable[sumTable==1] = "prev"
  sumTable <- as.matrix(sumTable)
  rownames(sumTable) <- paste0("Sequence", seq(nrow(sumTable)))

  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
                gp = gpar(fill = "#CCCCCC", col = NA))
    },
    # big blue
    prev = function(x, y, w, h) {
      grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
                gp = gpar(fill = "#DB4B71", col = NA))
    }
  )

  column_title = "OncoPrint"
  oncoPrint(sumTable[1:40, ],
           alter_fun = alter_fun, col = "#DB4B71",
           column_title = column_title)
  return(subTable3)
}
