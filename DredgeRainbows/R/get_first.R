#' Compute vector with occasions of first capture
#'
#' Copied from Kery and Schaub, Bayesian Population Analysis, 2012, pg 273
#'
#' @param x A capture history matrix
#'
#' @return A vector
#'
#' @export
get_first <- function(x) min(which(x!=0))