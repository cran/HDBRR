require(numDeriv)
require(parallel)
require(bigstatsr)
require(MASS)

.sysinfo <- Sys.info()
.OS <- .sysinfo[['sysname']]
.mc1.cores = getOption("mc.cores", 2L)
if(.OS == "Windows"){
  .mc1.cores = 1L
}

.mcsapply<-function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, mc.preschedule = TRUE,
                    mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = .mc1.cores,
                    mc.cleanup = TRUE, mc.allow.recursive = TRUE, affinity.list = NULL )
{
  answer <- mclapply(X = X, FUN = FUN, ...,mc.preschedule = mc.preschedule,
                     mc.set.seed = mc.set.seed, mc.silent = mc.silent, mc.cores = mc.cores,
                     mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive, affinity.list = affinity.list)
  if (USE.NAMES && is.character(X) && is.null(names(answer)))
    names(answer) <- X
  if (!isFALSE(simplify) && length(answer))
    simplify2array(answer, higher = (simplify == "array"))
  else answer
}

.big_fat_svd <- function(x,xxt=NULL){

  # Compute t(x) * x
  if(is.null(xxt)){
    x <- big_transpose(x)
    xxt <- big_crossprodSelf(x, big_scale(center = FALSE, scale = FALSE))
  }
  # Get v and d where a = u * d * t(v) the SVD of a
  v <- FBM(ncol(x), ncol(x))
  eig <- eigen(xxt[]+diag(1E-6,ncol(x)))
  v[] <- eig$vectors
  d <- sqrt(eig$values)

  # Get u, it will be of the same size of x
  # so that I store it as a FBM.
  u <- FBM(nrow(x), ncol(x))
  big_apply(u, a.FUN = function(X, ind, x, v, d) {
    X[ind, ] <- sweep(x[ind, ] %*% v[], 2, d, "/")
  }, a.combine = 'c', ind = rows_along(u),
  x = x, v = v, d = d)
  return(list(d=d,u=v,v=u,xxt=xxt))
}


.big_fat_tqr <- function(x,qmat=TRUE,xxt=NULL){
  x <- big_transpose(x)

  # Compute t(x) * x
  if(is.null(xxt)){
    xxt <- big_crossprodSelf(x, big_scale(center = FALSE, scale = FALSE))
  }

  # Get q and r where a = q*r the qr of a
  r <- FBM(ncol(x), ncol(x))
  rinv <- FBM(ncol(x), ncol(x))
  r[] <- chol(xxt[])
  if(qmat){
    rinv[] <- backsolve(r[],diag(ncol(x)))
    # Get q, it will be of the same size of x
    # so that I store it as a FBM.
    q <- FBM(nrow(x), ncol(x))
    big_apply(q, a.FUN = function(X, ind, x, rinv) {
      X[ind, ] <- sweep(x[ind, ] %*% rinv[],2,1,"*")
    }, a.combine = 'c', ind = rows_along(q),
    x = x, rinv = rinv)
  }
  return(list(q=q,r=r,xxt=xxt))
}


matop <- function(y,X,method = c("svd","qr"), bigmat = TRUE){
  n <- nrow(X)
  p <- ncol(X)
  Xfbm <- FBM(n, p)
  Xfbm[] <- X
  if(p > n && bigmat){
  	if(method[1]!="svd")
  	{
  	qr.X <- .big_fat_tqr(Xfbm,qmat=FALSE)
  	#q <- qr.X$q[]
  	r <- qr.X$r[]
  	}else{
  svd.X <- .big_fat_svd(Xfbm)
  ev <- svd.X$d^2
  D <-  svd.X$d
  L <- svd.X$u
  R <- svd.X$v
  }
  }

  if(p <= n || !bigmat){
    if(method[1] != "svd")
    {
      qr.X <- qr(X)
      r <- qr.X$qr
    }else{
      svd.X <- svd(X)
      ev <- svd.X$d^2
      D <-  svd.X$d
      L <- svd.X$u
      R <- svd.X$v
      Lfbm <- FBM(nrow(L),ncol(L))
      Rfbm <- FBM(nrow(R),ncol(R))
      Lfbm[] <- L
      Rfbm[]  <- R
      L <- Lfbm
      R <- Rfbm
    }
  }

 if(method[1]=="svd"){
  Ly <- big_cprodMat(L,as.matrix(y,ncol=1))
  return(list(y=y,X=X,D=D,L=L,R=R,ev=ev,Ly=Ly[],n=n,p=p))
  } else{
  	return(list(y=y,X=X,R=r,n=n,p=p))
  	}
}
############################################################


.logdn0 <- function(u,a1,a2,b1,b2,svdx){
  y <- svdx$y
  ev <- svdx$ev
  D <-  svdx$D
  L <- svdx$L
  R <- svdx$R
  Ly <- svdx$Ly
  n <- svdx$n
  p <- svdx$p
  cte <- -n/2*log(2*pi)+(a1/2)*log(b1/2)+(a2/2)*log(b2/2)-lgamma(a1/2)-lgamma(a2/2)+lgamma((n+a1+a2)/2)+(n+a1+a2)/2*log(2)
  ldn <- cte-1/2*sum(log((1-u)*ev+u))+(n-p)/2*u*(p > n)
  ldn <- ldn+((a1+n)/2-1)*log(1-u)
  if(n > p){
    W <- diag(D/(ev+u/(1-u)))
    WLy <- crossprod(W,Ly)
    betahat <- crossprod(t(R[]),WLy[])
    yhat <- crossprod(t(L[]),D*WLy[])
    SSE <- sum((y-yhat)^2)
    fu <- (1-u)*(SSE+b1)+u*(as.numeric(sum(betahat^2))+b2)
    ldn <- ldn+((a2+p)/2-1)*log(u)
  }
  else {
    fu <- u*(1-u)*sum(Ly^2/(ev*(1-u)+u))+(1-u)*b1+u*b2
    ldn <- ldn+((a2+n)/2-1)*log(u)
  }
  ldn <- ldn-(n+a1+a2)/2*log(fu)
  return(ldn)
}



.logdn <- function(u,a1,a2,b1,b2,svdx){
  y <- svdx$y
  ev <- svdx$ev
  D <-  svdx$D
  L <- svdx$L
  R <- svdx$R
  Ly <- svdx$Ly
  n <- svdx$n
  p <- svdx$p
  W <- diag(D/(ev+u/(1-u)))
  WLy <- tcrossprod(W,t(Ly))#big_cprodMat(W,Ly)
  yhat <- tcrossprod(L,t(D*WLy))#big_prodMat(L,D*WLy[])
   R <- R[,1:min(n,p)]
 betahat <- tcrossprod(R,t(WLy))#big_prodMat(Rfbm,WLy[])
  if(p>n){
    ev2 <- c(D^2+u/(1-u),rep(u/(1-u),p-n))
    fu <- u*(1-u)*sum(Ly^2/(ev*(1-u)+u))+(1-u)*b1+u*b2
  } else {
    ev2 <- D^2+u/(1-u)
    SSE <- sum((y-yhat)^2)
    fu <- (1-u)*(SSE+b1)+u*(as.numeric(sum(betahat^2))+b2)
  }
  invxx <- sapply(1:p,function(i)sum(R[i,]^2/ev2[1:min(n,p)]))
  vhat <- (n+a1+a2-1)/fu
  varb <- fu/((n+a1+a2-2)*(1-u))*invxx
  sigsqhat <- 1/(vhat*(1-u))
  sigbsqhat <- 1/(vhat*u)
  return(list(betahat=betahat,yhat=yhat,varb=varb,sigsqhat=sigsqhat,sigbsqhat=sigbsqhat,vhat=vhat,invxx=invxx))
}



.postu <- function(u,a1=1E-4,a2=1E-4,b1=1E-4,b2=1E-4,svdx, ncores){
  s <- length(u)
  n <- svdx$n
  p <- svdx$p
  dn <- rep(0,s)
  betahat <- matrix(0,nrow=s,ncol=p)
  varb <- matrix(0,nrow=s,ncol=p)
  invxx <- matrix(0,nrow=s,ncol=p)
  yhat <- matrix(0,nrow=n,ncol=s)
  sigsqhat <- rep(0,s)
  sigbsqhat <- rep(0,s)
  vhat <- rep(0,s)


  # Initiate cluster
  #cl <- makeCluster(ncores)
  #clusterExport(cl, ".logdn")
  #z <- seq(1:s)
  #lpostu <- parSapply(cl,z,function(i){
  #	.logdn(u[i],a1,a2,b1,b2,svdx)})
  #  stopCluster(cl)

  lpostu <- .mcsapply(1:s,function(i){
    .logdn(u[i],a1,a2,b1,b2,svdx)},mc.cores=ncores)

    betahat <- matrix(unlist(lpostu["betahat",]),ncol=p,byrow=T)
    yhat <- matrix(unlist(lpostu["yhat",]),ncol=n,byrow=T)
    varb <- matrix(unlist(lpostu["varb",]),ncol=p,byrow=T)
    sigsqhat <- matrix(unlist(lpostu["sigsqhat",]),ncol=1)
    sigbsqhat <- matrix(unlist(lpostu["sigbsqhat",]),ncol=1)
    vhat <- matrix(unlist(lpostu["vhat",]),ncol=1)
    invxx <- matrix(unlist(lpostu["invxx",]),ncol=p,byrow=T)

return(list(betahat=betahat,yhat=yhat,varb=varb,sigsqhat=sigsqhat,sigbsqhat=sigbsqhat,vhat=vhat,invxx=invxx))
}


HDBRR <- function(y,X,a1 = 1,a2 = 1,b1 = 1,b2 = 1/ncol(X),intercept = TRUE,npts = NULL,c = NULL,
                  corpred = NULL ,method = c("svd","qr"),bigmat = TRUE, ncores = 2){
  sysinfo <- Sys.info()
  OS <- sysinfo[['sysname']]
  if(OS == 'Windows'){
    ncores = 1
  }
  if(is.null(npts)) npts <- 200
  n <- nrow(X)
  p <- ncol(X)
  if(intercept) X <- cbind(rep(1,n),X)
  if(!is.vector(y)) y <- as.vector(y)
  Xall <- X
  yall <- y
  isNA <- is.na(y)
  whichNa <- which(isNA)
  y <- y[!isNA]
  X <- X[!isNA,]
  svdx <- matop(y,X,bigmat = bigmat)
  u0 <- seq(1E-6,1-1E-6,length=10)
  dn0 <- sapply(1:10,function(i).logdn0(u0[i],a1,a2,b1,b2,svdx))
  u0 <- u0[which.max(dn0)]
  par0 <- suppressWarnings(optim(u0,.logdn0,method="L-BFGS-B",lower=1E-6,upper=1-1E-10,a1=a1,a2=a2,b1=b1,b2=b2,svdx=svdx,control=list(fnscale=-1),hessian=F))
  umode <- par0$par
  hess <- suppressWarnings(hessian(.logdn0,umode,a1=a1,a2=a2,b1=b1,b2=b2,svdx=svdx))

  if(!is.nan(hess)){
    sd0 <- 1/sqrt(-hess)
    a <- umode*(umode*(1-umode)/sd0^2-1)
    b <- umode*(1-umode)/sd0^2-a-1
    u <- seq(qbeta(0.001,a,b),qbeta(0.999,a,b),length=npts)
  }else{
  	 u0 <- seq(max(0.001,umode-1/2),min(umode+1/2,0.999),length = npts)
     dn <- sapply(1:npts,function(i).logdn0(u0[i],a1,a2,b1,b2,svdx))
     idint <- range(which(exp(dn-max(dn))>0.001))
     u <- seq(u0[idint[1]],u0[idint[2]],length = npts)
   }

  dn <- sapply(1:npts,function(i).logdn0(u[i],a1,a2,b1,b2,svdx))
  dn <- exp(dn-max(dn))
  margu <- .postu(u,a1,a2,b1,b2,svdx,ncores = ncores)
  step <- 1/npts
  probs  <- step*dn
  probs <- probs/sum(probs)
  sigsqb <- margu$varb
  betahat <- c(probs%*%margu$betahat)
  varb <-probs%*%margu$varb
  varb <- c(t(varb)+(t(margu$betahat)-c(betahat))^2%*%probs)
  yhat <- if(sum(isNA)>0) c(Xall%*%betahat) else c(probs%*%margu$yhat)#big_prodMat(Xall,matrix(betahat,ncol=1))
  sigsqhat <- c(probs%*%margu$sigsqhat)
  sigbsqhat <- c(probs%*%margu$sigbsqhat)
  uhat <- c(u%*%probs)


if(!is.null(c)){
  vhat <- 1/sigsqhat+1/sigbsqhat
  lambda <- c(uhat/(1-uhat))
  SSE <- sum((y-yhat[!isNA])^2)
  ssx <- sapply(1:(p+intercept),function(i)sum(X[,i]^2))
  beta.cond <- sapply(1:(p+intercept),function(i)t(X[,i])%*%(y-yhat[!isNA]+X[,i]*betahat[i])/(ssx[i]+lambda))
  phi <- 1/2
  phat <- 1-(1-phi)*sqrt(lambda/(ssx+lambda))*exp(-vhat*(1-uhat)*ssx*betahat^2/2)*exp(-uhat*vhat*betahat^2/2)/exp(-vhat*(1-uhat)*(ssx+uhat/(1-uhat))*(betahat-beta.cond)^2/2)
  delta <- sqrt(sigbsqhat*2*log(c)*c^2/(c^2-1))
 } else {
 	  phat <- NULL
 	  delta <- NULL
 	  }

if(method[1] == "svd" && !is.null(corpred)){
	if(corpred != "eb"){
      num1 <- sapply(1:length(u),function(i)(1-u[i])*svdx$ev)
      pj <- num1/(num1+u)
      wgt <- probs
  }
  else{
 	  num1 <- (1-uhat)*svdx$ev
 	  pj <- num1/(num1+uhat)
 	  wgt <- 1
 	}
  usqrt <- svdx$L[]^2
  num2 <- usqrt%*%pj
  den2 <- usqrt%*%pj^2
  corr <- apply(num2/sqrt(den2)*wgt,1,sum)
  edf <- sum(pj*wgt)
}
else{
 	corr <- edf <- NULL
 }
 return(list(betahat=betahat,yhat=yhat,varb=varb,sigsqhat=sigsqhat,sigbsqhat=sigbsqhat,u=u,postu=probs,uhat=uhat,umode=umode,whichNa=whichNa,phat=phat,delta=delta,edf=edf,corr=corr))
}




