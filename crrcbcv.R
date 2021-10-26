#' Bias-Corrected Variance for Competing Risks Regression with Clustered Data
#'
#' @description
#' Small-sample bias-corrected variance for regression modeling using proportional subdistribution hazards with
#' clustered right censored data. (Zhou et al., 2012) Failure times within the same cluster are dependent.
#'
#' Four types of bias correction are included: the MD-type correction by Mancl and DeRouen (2001),
#' the KC-type correction by Kauermann and Carroll (2001), the FG-type correction by Fay and Graubard (2001),
#' and the MBN-type correction by Morel, Bokossa, and Neerchal (2003).
#'
#' @import abind crrSC pracma survival stats
#'
#' @param beta the estimated regression coefficients from `crrc`
#' @param ftime vector of failure/censoring times
#' @param fstatus vector with a unique code for each failure type and a separate code for censored observations
#' @param cov1 matrix (nobs x ncovs) of fixed covariates (either cov1, cov2, or both are required)
#' @param cov2 matrix of covariates that will be multiplied by functions of time;
#' if used, often these covariates would also appear in cov1 to give a prop hazards effect plus a time interaction
#' @param tf functions of time. A function that takes a vector of times as an argument and returns a matrix whose
#' jth column is the value of the time function corresponding to the jth column of cov2 evaluated at the
#' input time vector. At time `tk`, the model includes the term `cov2[,j]*tf(tk)[,j]` as a covariate
#' @param cluster clustering indicator
#' @param failcode code of fstatus that denotes the failure type of interest
#' @param cencode code of fstatus that denotes censored observations
#' @param subset a logical vector specifying a subset of cases to include in the analysis
#' @param na.action a function specifying the action to take for any cases missing any of
#' ftime, fstatus, cov1, cov2, cengroup, or subset
#' @param var.type a string or a vector of strings with value(s) selected from \{`"MD"`, `"KC"`, `"FG"`, `"MBN"`\}
#'
#' @return
#' Returns a list of class crr, with components corresponding to `var.type`
#' \describe{
#'   \item{`$MD`}{the MD-type bias-corrected variance covariance matrix for `beta`}
#'   \item{`$KC`}{the KC-type bias-corrected variance covariance matrix for `beta`}
#'   \item{`$FG`}{the FG-type bias-corrected variance covariance matrix for `beta`}
#'   \item{`$MBN`}{the MBN-type bias-corrected variance covariance matrix for `beta`}
#' }
#'
#' @export
#'
#' @author
#' Xinyuan Chen, `<xchen@math.msstate.edu>`
#'
#' Fan Li, `<fan.f.li@yale.edu>`
#'
#' @references
#' Chen X, Li F. (2022). Finite-sample adjustments in variance estimators for
#' clustered competing risks regression. Statistics in Medicine. 00(Under Review): 1-24.
#'
#' Zhou B, Fine J, Latouche A, Labopin M. (2012). Competing risks regression for clustered Data. Biostatistics.
#' 13(3): 371-383.
#'
#' @seealso `crrSC`
#'
#' @examples
#' library(crrcbcv)
#' data(cls)
#' mod.est = crrc(ftime=cls$T_obs, fstatus=cls$eps, cov1=cls[,c('X_1','X_2')], cluster=cls$I)
#' crrcbcv(beta=mod.est$coef, ftime=cls$T_obs, fstatus=cls$eps, cov1=cls[,c('X_1','X_2')],
#' cluster=cls$I, var.type=c('MD','KC','FG','MBN'))
crrcbcv = function(beta, ftime, fstatus, cov1, cov2, tf, cluster, failcode=1,
                   cencode=0, subset, na.action=na.omit, var.type="MD") {

  d = data.frame(ftime=ftime,fstatus=fstatus,
                 cluster=if (missing(cluster)) 1:length(fstatus) else cluster)

  if (!missing(subset)) {
    d = d[subset,]
    if (!missing(cov1)) {cov1 = cov1[subset,]}
    if (!missing(cov2)) {cov2 = cov2[subset,]}
  }

  tmp = nrow(d)
  d = na.action(d)

  if (nrow(d) != tmp) {cat(format(tmp-nrow(d)),'cases omitted due to missing values\n')}

  if (!missing(cov1)) {
    cov1 = as.matrix(cov1)
    np1 = ncol(cov1)
    cov1 = na.action(cov1)
  } else {np1 = 0}

  if (!missing(cov2)) {
    cov2 = as.matrix(cov2)
    np2 = ncol(cov2)
    cov2 = na.action(cov2)
  } else {np2 = 0}

  np = np1 + np2

  if (np == 0) {return("No covariates included.")} # breaks out of the function

  # sort observations by time
  Y = d$ftime
  b = order(Y)
  Y = sort(Y)
  ID = d$cluster
  ID = ID[b]
  d = d[b,]

  UID = sort(unique(ID))
  IDind = pracma::zeros(length(UID), length(Y))
  for (i in 1:length(UID)){IDind[i, ID==UID[i]] = 1}
  n = length(UID) # number of clusters
  ny = length(Y)

  cenind = ifelse(d$fstatus==cencode,1,0)
  fstatus = ifelse(d$fstatus==failcode,1,0)

  # increment of time
  dY = Y-c(0,Y[1:(ny-1)])
  # trick to obtain at-risk process
  # each row is an individual
  # each column is a specific time point
  IndYY = (t(pracma::repmat(Y,ny,1)) >= pracma::repmat(Y,ny,1))
  IndYY_oc = IndYY
  oc = !(d$fstatus %in% c(cencode, failcode))
  IndYY_oc[oc,] = 1

  if (np2 == 0) {
    cov1 = as.matrix(cov1)
    cov1 = cov1[b,]
    cov2 = 0
    # X is an array, the first dimension is person,
    # the second dimension is time, the third dimension is covariate dimensionality
    X = abind::abind(lapply(1:ncol(cov1), function(x) cov1[,x] %o% rep(1,ny)), along=3)
  } else if (np1 == 0) {
    cov2 = as.matrix(cov2)
    cov2 = cov2[b,]
    cov1 = 0
    X = abind::abind(lapply(1:ncol(cov2), function(x) cov2[,x] %o% tf(Y)), along=3)
  } else {
    cov1 = as.matrix(cov1)
    cov1 = cov1[b,]
    cov2 = as.matrix(cov2)
    cov2 = cov2[b,]
    X1 = abind::abind(lapply(1:ncol(cov1), function(x) cov1[,x] %o% rep(1,ny)), along=3)
    X2 = abind::abind(lapply(1:ncol(cov2), function(x) cov2[,x] %o% tf(Y)), along=3)
    X = abind::abind(X1, X2, along=3)
  }

  # the rate of counting process of event, or dN(t)
  NN = diag(fstatus)

  # estimate IPCW
  u = survival::survfit(survival::Surv(Y,cenind)~1)
  u = summary(u, times=sort(Y*(1-.Machine$double.eps)))
  uuu = u$surv

  # a matrix of weights at each t
  # each row is an individual
  # each column is a time
  Weights_tmp = pracma::repmat(uuu,length(uuu),1) / sapply(uuu, function(x) x*(x>=uuu)+uuu*(x<uuu))
  res = pracma::repmat(rep(1,length(uuu)), length(uuu), 1)
  for (k1 in which(cenind==1)) {res[k1,-(1:k1)] = 0}
  Weights = res*Weights_tmp

  # Compute dHc; this requires the real IndYY, not IndYY_oc
  NNc = diag(cenind)
  dHc = c(-log(uuu)[1], diff(-log(uuu)))

  Xbeta = matrix(0,ny,ny)
  for(k in 1:np) {Xbeta = Xbeta + X[,,k]*beta[k]}

  # S0, S1, S2
  S0 = n**(-1) * colSums(Weights*IndYY_oc*exp(Xbeta))
  S1 = matrix(0, nrow=np, ncol=ny)
  for(k in 1:np){
    S1[k,] = n**(-1) * colSums(Weights*IndYY_oc*X[,,k]*exp(Xbeta))
  }
  S2 = array(0, c(ny,np,np))
  for(k in 1:np){
    for(s in 1:np){
      S2[,k,s] = n**(-1) * colSums(Weights*IndYY_oc*X[,,k]*X[,,s]*exp(Xbeta))
    }
  }

  # Estimation of marginal and conditional baseline hazard
  # dHY is baseline hazard function (not individual martigale)
  # HY is cumulative baseline hazard function (not individual martigale)
  dHY = n**(-1) * colSums(NN*Weights / pracma::repmat(S0,ny,1))
  HY = cumsum(dHY)

  # dM_hat # individual martingale increment
  dMij = NN - IndYY_oc*exp(Xbeta) * pracma::repmat(dHY,ny,1)

  # compute dMic
  dMijc = NNc - IndYY * pracma::repmat(dHc,ny,1)
  XmbarX = array(0, c(ny,ny,np))
  for (k in 1:np) {
    XmbarX[,,k] = X[,,k] - pracma::repmat(S1[k,]/S0,ny,1)
  }

  # compute eta and psi; eta is easier
  eta_ij = matrix(0,ny,np)
  for (k in 1:np) {
    eta_ij[,k] = rowSums((X[,,k] - pracma::repmat(S1[k,]/S0,ny,1)) * Weights*dMij)
  }

  # compute psi is is more difficult
  # compute q(u)
  qu = matrix(0,nrow=ny,ncol=np)
  for (i in 1:ny) {
    gr = sum(Y[i] > Y)
    ln = sum(Y[i] <= Y)
    if (gr == 0) {
      qu[i,] = rep(0, np)
    } else {
      q_ind = cbind(matrix(0,ncol=gr,nrow=gr), matrix(1,ncol=ln,nrow=gr))
      qtemp1 = q_ind*Weights[1:gr,]*dMij[1:gr,]
      qu[i,] = -n**(-1) * do.call(c, lapply(1:np, function(x) sum(XmbarX[1:gr,,x]*qtemp1)))
    }
  }

  # compute pi(u)
  pi_u = n**(-1) * colSums(IndYY)

  # compute psi_ij
  psi_ij = do.call(cbind, lapply(1:np, function(x) rowSums(pracma::repmat((qu/pi_u)[,x],ny,1) * dMijc))) # ny*np

  eta_i = as.matrix(stats::aggregate(eta_ij, by=list(ID), FUN=sum)[,-1])
  psi_i = as.matrix(stats::aggregate(psi_ij, by=list(ID), FUN=sum)[,-1])

  U_i = eta_i + psi_i

  Omega_ij = array(0, c(ny,np,np))
  for (k in 1:np) {
    for (s in 1:np) {
      Omega_ij[,k,s] = rowSums(pracma::repmat(S2[,k,s]/S0-(S1[k,]*S1[s,])/(S0**2),ny,1) * Weights*dMij) +
        rowSums((X[,,k]-pracma::repmat(S1[k,]/S0,ny,1))*X[,,s]*IndYY_oc*Weights*exp(Xbeta) *
                  pracma::repmat(dHY,ny,1))
    }
  }

  Omega_i = array(0,c(n,np,np))
  for(i in 1:length(UID)){
    for (k in 1:np){
      for (s in 1:np){
        Omega_i[i,k,s] = sum(Omega_ij[IDind[i,]==1,k,s])
      }
    }
  }

  Vm = solve(apply(Omega_i, c(2,3), sum))
  Vs = Vm %*% crossprod(U_i) %*% Vm

  nomMD = matrix(0,n,np)
  H_iU_i = matrix(0,n,np)
  for (i in 1:n) {
    nomMD[i,] = solve(diag(np)-Omega_i[i,,]%*%Vm) %*% U_i[i,]
    H_iU_i[i,] = diag(1/sqrt(1-pmin(0.75,c(diag(Omega_i[i,,] %*% Vm))))) %*% U_i[i,]
  }

  outlist = list()

  # compute MD type correction variance
  if ("MD" %in% var.type) {
    BC_MD = Vm %*% crossprod(nomMD) %*% Vm
    outlist[["MD"]] = BC_MD
  }

  # compute KC type correction variance
  if ("KC" %in% var.type) {
    BC_KC = 0.5*(Vm %*% crossprod(nomMD,U_i) %*% Vm + Vm %*% crossprod(U_i,nomMD) %*% Vm)
    outlist[["KC"]] = abs(diag(BC_KC)) + (matrix(1,ncol=np,nrow=np)-diag(np))*BC_KC
  }

  # compute FG type correction variance
  if ("FG" %in% var.type) {
    BC_FG =  Vm %*% crossprod(H_iU_i) %*% Vm
    outlist[["FG"]] = BC_FG
  }

  # compute MBN type correction variance
  if ("MBN" %in% var.type) {
    BC_MBN = (ny-1)*n/((ny-np)*(n-1))*Vs +
      min(0.5,np/(n-np))*max(1,sum(diag((ny-1)*n/((ny-np)*(n-1))*(Vm %*% crossprod(U_i))))/np)*Vm
    outlist[["MBN"]] = BC_MBN
  }

  class(outlist) = 'crrcbcv'
  return(outlist)
}

