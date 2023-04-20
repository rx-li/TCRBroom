#' Identify outliers 
#'
#' @description The function is used to identify the outliers
#'
#' @param vec   Named vector of the pairwise   
#' @param iters Number of iterations
#' @param c     Quantile threshold
#' @param jags_model The hierarchical model 
#' @param jags_params Parameters of the hierarchical model 
#' @param jags_inits Inits of the parameters
#'
#' @import R2jags 
#'
#'
#' @export
#'
#'
# data: data used for jags
detectOutliers <- function(vec, iters, c, jags_inits, jags_params, jags_model){
  quantile = c
  q.l <- NULL
  index.l <- NULL
  coe.l <- NULL
  k.l <- NULL
  
  y <- vec 
  k <- length(y)
  
  i <- 0
  for(i in seq(iters)){
    jagsfit <- jags(data = list("y", "k"),
                    inits = jags_inits, jags_params,
                    n.iter = 5000, model.file = jags_model)
    jagsfit_coe <- jagsfit$BUGSoutput$mean
    q <- qbeta(quantile, jagsfit_coe$alpha, jagsfit_coe$beta)
    y <- vec[vec < q]
    
    if(k <= length(y)){
      cat("Converge! Iters =", i-1)
      break
    }
    
    k <- length(y)
    
    q.l <- c(q.l, q)
    index.l <- c(index.l, list(vec[vec >= q]))
    coe.l <- c(coe.l, list(jagsfit_coe))
    k.l <- c(k.l, k)
    i <- i + 1
  }
  return(list(quantiles = q.l, 
              outliers = index.l, 
              coefficients = coe.l))
}
