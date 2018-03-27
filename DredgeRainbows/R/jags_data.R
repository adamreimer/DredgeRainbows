#' Function to create dataset for Jags
#'
#' Copied from Kery and Schaub, Bayesian Population Analysis, 2012, pg 273
#'
#' @param PSI.STATE state process matrix
#' @param PSI.OBS Observation process matrix
#' @param marked number of marked individuals per event
#'
#' @return A list
#'
#' @export
jags_data <- function(PSI.STATE, PSI.OBS, marked){
  # Execute simulation function
  sim <- simul_ms(PSI.STATE, PSI.OBS, marked)
  CH <- sim$CH
  f <- apply(CH, 1, get_first)
  
  # Recode CH matrix: note, a 0 is not allowed in WinBUGS!
  # 1 = seen alive in A, 2 = seen alive in B, 3 = not seen
  rCH <- CH          # Recoded CH
  rCH[rCH==0] <- 3
  
  # Bundle data
  list(y = rCH,
       f = f,
       n.occasions = dim(rCH)[2],
       nind = dim(rCH)[1],
       z = known_state_ms(rCH, 3))
}