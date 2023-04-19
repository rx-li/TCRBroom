#' Identify crosscohort overlap outliers 
#'
#' @description The function is used to identify crosscohort overlap outliers.
#'
#' @param ref Reference vector
#' @param vec   Your own vector    
#'
#'
#' @export
#'
#'
# data: data used for jags
Pcalculator <- function(ref, vec){
  p <-  unlist(lapply(vec, function(x) sum(ref > x) / length(ref)))
  fdr <- p.adjust(p, method="BH")
  out <- data.frame(
    value=vec, p=p, p_adjust=fdr 
  )
}