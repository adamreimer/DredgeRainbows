#' Function to create initial values for unknown z
#'
#' Copied from Kery and Schaub, Bayesian Population Analysis, 2012, pg 276
#'
#' @param ch capture history matrix
#' @param f occassion of first capture
#'
#' @return A matrix
#'
#' @export
ms_init_z <- function(ch, f){
  for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
  states <- max(ch, na.rm = TRUE)
  known.states <- 1:(states-1)
  v <- which(ch==states)
  ch[-v] <- NA
  ch[v] <- sample(known.states, length(v), replace = TRUE)
  return(ch)
}