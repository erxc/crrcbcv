#' Clustered competing risks simulated data
#'
#' @name cls
#'
#' @docType data
#'
#' @description sample of 20 clusters with an average cluster size of 20
#'
#' @usage data(cls)
#'
#' @format
#' A data frame containing 20 clusters with an average cluster size of 20 and the following five variables.
#' Simulation is detailed in the paper Finite-sample adjustments in variance estimators for
#' clustered competing risks regression. Chen, Li. 2022. Under Review. Statistics in Medicine.
#' \describe{
#'   \item{`I`}{id of clusters}
#'   \item{`X_1`}{a cluster-level covariate generated from the standard normal distribution}
#'   \item{`X_2`}{an individual-level covariate generated from the standard normal distribution}
#'   \item{`eps`}{event type. 0=censored, 1, 2}
#'   \item{`T_obs`}{observed event time}
#' }
#'
#' @keywords data
#'
#' @author
#' Xinyuan Chen, `<xchen@math.msstate.edu>`
#'
#' Fan Li, `<fan.f.li@yale.edu>`
#'
#' @references Chen X, Li F. (2022). Finite-sample adjustments in variance estimators for clustered competing risks
#' regression. Statistics in Medicine. 00(Under Review): 1-24.
#'
#' @examples
#' data(cls)
NULL
