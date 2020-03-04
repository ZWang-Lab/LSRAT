Get.p.eigen<-function(Q, A1){
  #here
  lambda <- Get_Lambda(A1);
  param<-Get_Liu_Params_Mod_Lambda(lambda); 
  # print(sum(param$muQ))
  # print(param$muQ)
  Q.Norm<-(Q - param$muQ)/param$sigmaQ;Q.Norm<-Q.Norm * param$sigmaX + param$muX;
  p.value<- pchisq(Q.Norm,  df = param$l,ncp=param$d, lower.tail=FALSE)
  return(p.value);
}

Get.p.eigen.burden<-function(Q,  V){
  # print(V)
  p.value  <- pchisq(Q/V, df = 1, lower.tail = FALSE)
  return(p.value);
}


Get.p<-function(Q,B.Q){
  B.mean<-mean(B.Q);
  # print(B.mean)
  B.variance<-var(B.Q)
  B.kurtosis <- mean((t(B.Q)-B.mean)^4)/B.variance^2-3
  B.df<-(B.kurtosis>0)*12/B.kurtosis+(B.kurtosis<=0)*10000000
  B.p<-t(1-pchisq((t(Q)-B.mean)*sqrt(2*B.df)/sqrt(B.variance)+B.df,B.df))
  return(B.p)
}


Get.p.eigen.skato2 <- function(score, V, rho, method = "saddle") {
    n.r <- length(rho);lambdas <- vector("list", n.r);pval <- qval <- rep(NA, n.r)
    n.p <- length(score); U  <- score;
    Q <- (1-rho)*sum(U^2)+rho*sum(U)^2
    Burden.score <- Burden.var <- Burden.pval <- SKAT.pval <- NA
    for(i in 1:n.r) {
    if(rho[i]==1) {
      Burden.score <- sum(U)
      Burden.var <- sum(V)
      Burden.pval <- pchisq(Burden.score^2/Burden.var, df=1, lower.tail=FALSE)
      lambdas[[i]] <- Burden.var
      pval[i] <- Burden.pval
      next
    }
  if(rho[i]!=0) {
      R.M <- matrix(rho[i], n.p, n.p)
      diag(R.M) <- 1
      R.M.chol <- t(chol(R.M, pivot = TRUE))
      V.temp <- crossprod(R.M.chol, crossprod(V, R.M.chol))
  } else V.temp <- V
   lambda <- eigen(V.temp, only.values = TRUE, symmetric=TRUE)$values
      lambdas[[i]] <- lambda[lambda > 0]
  pval[i] <- .Q_pval(Q[i], lambdas[[i]], method = method)
  if(rho[i]==0) SKAT.pval <- pval[i]
    }
    minp <- min(pval)
    for(i in 1:n.r) {
  df <- sum(lambdas[[i]]^2)^2/sum(lambdas[[i]]^4)
  qval[i] <- (qchisq(minp, df, lower.tail = FALSE)-df)/sqrt(2*df)*sqrt(2*sum(lambdas[[i]]^2))+sum(lambdas[[i]])
    }
    ZMZ <- tcrossprod(rowSums(V))/sum(V)
    V.temp <- V - ZMZ
    lambda <- eigen(V.temp, only.values = TRUE, symmetric = TRUE)$values
    lambda <- lambda[lambda > 0]
    muq <- sum(lambda)
    varq <- sum(lambda^2) * 2 + sum(ZMZ * V.temp) * 4
    df <- sum(lambda^2)^2/sum(lambda^4)
    tau <- rho * sum(V) + sum(V %*% V)/sum(V) * (1 - rho)
    re <- tryCatch({
        integrate(function(x){
          t1 <- tau %x% t(x)
          re<-pchisq((apply((qval - t1)/(1-rho),2,min) - muq)/sqrt(varq)*sqrt(2*df) + df, df=df) * dchisq(x,df=1)
          return(re)
  }, lower = 0, upper = 40, subdivisions = 2000, abs.tol = 10^-25)
    }, error=function(e) NA)
    return(list(p = min(1-re[[1]], minp*n.r), minp = minp, minp.rho = rho[which.min(pval)],
    Burden.score=Burden.score, Burden.var=Burden.var, Burden.pval=Burden.pval,
    SKAT.pval=SKAT.pval))
}


 Get.p.eigen.skato <- function(Q1,Q2, V, rho, burden.pval) {
  n.r <- length(rho);lambdas <- vector("list", n.r);pval <- rep(NA, n.r)
  n.p <- ncol(V);U = Q1*rho + Q2*(1-rho);
  for(i in 1:n.r){
    if(rho[i]==0){
      pval[i]=burden.pval
    }else{
      R.M <- matrix(1-rho[i], n.p, n.p);diag(R.M) <- 1
      R.V = V%*%R.M
      pval[i]  <- .quad_pval(U[i], R.V)
    }
  }
  minp = min(pval)
  min.rho = rho[which.min(pval)]
  return(list(minp = minp, min.rho = min.rho))
 }


#utility function from GMMAT

.Q_pval <- function(Q, lambda, method = "davies") {
    if(method == "davies") {
        tmp <- suppressWarnings(CompQuadForm::davies(q = Q, lambda = lambda, acc = 1e-6))
        pval <- tmp$Qq
        if((tmp$ifault > 0) | (pval <= 1e-5) | (pval >= 1)) method <- "liu"
    }

    if(method == "kuonen") {
      pval <- .pKuonen(x = Q, lambda = lambda)
      if(is.na(pval)) method <- "liu"
    }

    if(method == "saddle") {
      pval <- saddle(x = Q, lambda = lambda)
      if(is.na(pval)) method <- "liu"
    }

    if(method == "liu") pval <- CompQuadForm::liu(q = Q, lambda = lambda)
      return(pval)

}

.quad_pval <- function(Q, V, method = "davies") {
    lambda <- eigen(V, only.values = TRUE, symmetric=TRUE)$values
    lambda <- lambda[lambda > 0]
    # print(sum(lambda))
    pval <- .Q_pval(Q, lambda, method = method)
    return(pval)
}


.smmat_pvalue <- function(U, V, p.burden){
      V.rowSums <- rowSums(V)
      U <- U - V.rowSums * sum(U) / sum(V)
      V <- V - tcrossprod(V.rowSums) / sum(V)
      Q = sum(U^2)
      if(mean(abs(V)) < sqrt(.Machine$double.eps)) p.smmat <- p.burden
    else p.smmat <- tryCatch(pchisq(-2*log(p.burden)-2*log(.quad_pval(Q = Q, V = V)), df = 4, lower.tail = FALSE), error = function(e) { p.burden })
    return(p.smmat)
}    



.smmat_pvalue_retro <- function(Q1.e, V.e, p.burden){
  
  p.smmat <- tryCatch(pchisq(-2*log(p.burden)-2*log(.quad_pval(Q1.e, V.e)), df = 4, lower.tail = FALSE), error = function(e) { p.burden })
    return(p.smmat)
}    



snp_impute<-function(snp.mat, impute="mean")
{
	snp.imp <- snp.mat;

	for(i in 1:NCOL(snp.mat) )
	{
		s.mat.i <- snp.mat[,i] ;
		s.miss <- which( is.na(s.mat.i) );

		if (length(s.miss)>0)
		{
			if(impute=="mean")
			{
				s.mat.i[s.miss] <- mean(s.mat.i, na.rm=T);
			}
			else
			{
				n.s0 <- length( which( s.mat.i == 0 ) );
				n.s1 <- length( which( s.mat.i == 1 ) );
				n.s2 <- length( which( s.mat.i == 2 ) );
				n.s  <- length(s.mat.i)

				r.miss<- runif( length(s.miss) );
				r.snp <- rep(2, length(s.miss));
				r.snp[r.miss <= n.s0/n.s ]<-0;
				r.snp[r.miss <= (n.s0 + n.s1)/n.s ]<-1;
				s.mat.i[s.miss] <- r.snp;
			}
		}

		if (mean(s.mat.i)/2>0.5) s.mat.i <- 2 - s.mat.i;

		snp.imp[,i] <- s.mat.i;
	}

	return(snp.imp);
}





#utility function of LGEWIS
Get.sqrt<-function(A){
  a.eig <- eigen(A,symmetric=TRUE)
  ID1<-which(a.eig$values > 0)
  if(length(ID1)== 0){stop("Error to obtain matrix square!")}
  a.sqrt <- a.eig$vectors[,ID1] %*% diag(sqrt(a.eig$values[ID1])) %*% t(a.eig$vectors[,ID1])
  return(a.sqrt)
}



Get.inv.sqrt<-function(A){
  a.eig <- eigen(A,symmetric=TRUE)
  ID1<-which(a.eig$values > 0)
  if(length(ID1)== 0){stop("Error to obtain matrix square!")}
  a.sqrt <- a.eig$vectors[,ID1] %*% diag(sqrt(a.eig$values[ID1])^-1,length(ID1)) %*% t(a.eig$vectors[,ID1])
  return(a.sqrt)
}

Get.inverse<-function(A){
  a.eig <- eigen(A,symmetric=TRUE)
  ID1<-which(a.eig$values > 0)
  if(length(ID1)== 0){stop("Error to obtain matrix inverse!")}
  a.inverse <- a.eig$vectors[,ID1] %*% diag(a.eig$values[ID1]^-1) %*% t(a.eig$vectors[,ID1])
  return(a.inverse)
}


calAR1 = function(n.sample, n.time, rho){
H <- abs(outer(1:n.time, 1:n.time, "-"))
V <- rho^H
V_list = lapply(seq_len(n.sample), function(X) V)
AR1 = as.matrix(bdiag(V_list))
return(AR1)
}
calB = function(n.sample, n.time){
  M10 <- matrix(0, n.sample*n.time, n.sample*n.time)
  for(i in 1:n.sample) M10[(i-1)*n.time+(1:n.time), (i-1)*n.time+(1:n.time)] <- 1
  return(M10)
  }

inv.logit=function(x){
	return(exp(x)/(1+exp(x)))
}

logit =function(x){
    return (log(x/(1-x)))
}

Get_Lambda<-function(K)
{
	out.s <- try(eigen(K,symmetric=TRUE, only.values = TRUE))
	#print(out.s$values)

	if (class(out.s)=="try-error") { show(K); browser(); }

	#out.s1<-eigen(K,symmetric=TRUE)
	#print(out.s1$values)

	lambda1<-out.s$values
	IDX1<-which(lambda1 >= 0)

	# eigenvalue bigger than sum(eigenvalues)/1000
	IDX2<-which(lambda1 > mean(lambda1[IDX1])/100000)

	if(length(IDX2) == 0){
		stop("No Eigenvalue is bigger than 0!!")
	}

	lambda<-lambda1[IDX2]
	return(lambda)
}


Get_Liu_Params<-function(c1){
  ## Helper function for getting the parameters for the null approximation
  muQ<-c1[1]
  sigmaQ<-sqrt(2 *c1[2])
  s1 = c1[3] / c1[2]^(3/2)
  s2 = c1[4] / c1[2]^2

  beta1<-sqrt(8)*s1
  beta2<-12*s2
  type1<-0

  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 *a^3 - a^2
    l = a^2 - 2*d
  } else {
    type1<-1
    a = 1/s1
    d = 0
    l = 1/s1^2
  }
  muX <-l+d
  sigmaX<-sqrt(2) *a

  re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  return(re)
}

Get_Liu_Params_Mod_Lambda<-function(lambda){
  ## Helper function for getting the parameters for the null approximation

  c1<-rep(0,4)
  for(i in 1:4){
    c1[i]<-sum(lambda^i)
  }

  muQ<-c1[1]
  sigmaQ<-sqrt(2 *c1[2])
  s1 = c1[3] / c1[2]^(3/2)
  s2 = c1[4] / c1[2]^2

  beta1<-sqrt(8)*s1
  beta2<-12*s2
  type1<-0

  #print(c(s1^2,s2))
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 *a^3 - a^2
    l = a^2 - 2*d
  } else {
    type1<-1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX <-l+d
  sigmaX<-sqrt(2) *a

  re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  return(re)
}

Get_Satterthwaite<-function(muQ, varQ){

    a1<-varQ/muQ/2
    a2<-muQ/a1

    re<-list(df=a2, a=a1)
    return(re)
}

covSelector= function(index, n.time){
    cov.s = c()
    for(i in index){
    cov.s = c(cov.s, ((i-1)*n.time+1):(i*n.time))
    }
    return(cov.s)
}



### saddle point approximation ###
saddle = function(x,lambda){
  d = max(lambda)
  lambda = lambda/d
  x = x/d
  k0 = function(zeta) -sum(log(1-2*zeta*lambda))/2
  kprime0 = function(zeta) sapply(zeta, function(zz) sum(lambda/(1-2*zz*lambda)))
  kpprime0 = function(zeta) 2*sum(lambda^2/(1-2*zeta*lambda)^2)
  n = length(lambda)
  if (any(lambda < 0)) {
    lmin = max(1/(2 * lambda[lambda < 0])) * 0.99999
  } else if (x>sum(lambda)){
    lmin = -0.01
  } else {
    lmin = -length(lambda)/(2*x)
  }
  lmax = min(1/(2*lambda[lambda>0]))*0.99999
  hatzeta = uniroot(function(zeta) kprime0(zeta) - x, lower = lmin, upper = lmax, tol = 1e-08)$root
  w = sign(hatzeta)*sqrt(2*(hatzeta*x-k0(hatzeta)))
  v = hatzeta*sqrt(kpprime0(hatzeta))
  if(abs(hatzeta)<1e-4){
    return(NA)
  } else{
    return( pnorm(w+log(v/w)/w, lower.tail=FALSE) )
  }
}


colSds<-function(mat, na.rm=T)
{
  r<-c();
  for(i in 1:dim(mat)[2])
    r <- c(r, sd(mat[,i], na.rm=na.rm));
  return(r);
}
require("data.table")
library(plyr)
SNPs_maps<-function(bimfile, rsID){
  bim<-fread(bimfile)
  bim.matched<-bim[which(bim$V2 %in% rsID),c(2, 1, 4)]
  colnames(bim.matched)= c('SNP.name', 'chr', 'bp') 
  bim.matched$SNP.name<-as.character(bim.matched$SNP.name)
  bim.matched <- data.frame(bim.matched)
  return(bim.matched)
}

read.grm  = function(filename){
  mat <- scan(filename, what = numeric())
  ncol <- (sqrt(8 * length(mat) + 1) - 1) / 2
  diag_idx <- cumsum(seq.int(ncol))
  split_idx <- cummax(sequence(seq.int(ncol)))
  split_idx[diag_idx] <- split_idx[diag_idx] - 1
  splitted_rows <- split(mat, f = split_idx)
  mat_full <- suppressWarnings(do.call(rbind, splitted_rows))
  mat_full[upper.tri(mat_full)] <- t(mat_full)[upper.tri(mat_full)]
  return(mat_full)
}
