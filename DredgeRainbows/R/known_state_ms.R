#' Function to create known latent states z
#'
#' Copied from Kery and Schaub, Bayesian Population Analysis, 2012, pg 276
#'
#' @param ms capture history matrix
#' @param notseen label for ‘not seen'
#'
#' @return A matrix
#'
#' @export
known_state_ms <- function(ms, notseen){
  # notseen: label for ‘not seen’
  state <- ms
  state[state==notseen] <- NA
  for (i in 1:dim(ms)[1]){
    m <- min(which(!is.na(state[i,])))
    state[i,m] <- NA
  }
  return(state)
}