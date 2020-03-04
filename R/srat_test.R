#' Calculate prospective and retrospective P-values for GEE or GLMM model
#'
#' This function tests a SNPs for a given SNP set for a given lbrat estimated null model.
#'
#' @param srat.est The output of function "lsrat_est()" or "rsmmat_est()"
#' @param G The genotype matrix, an m*q matrix where m is the number of subjects and q is the total number genetic variants. 
#' @param impute.method choose the iputation method when there is missing genotype. Optional options are: 'random', 'fixed' or 'bestguess'.
#' @param GRM takes m-by-m genetic correlation matrix or kinship matrix.
#' @param tests a character vector indicating which LSRAT/RSMMAT tests should be performed ("B" for the burden test, "S" for SKAT, "C" for variant-level ACAT test, "O" for SKAT-O, "E" for the SMMAT test, "A" for omnibus ACAT test). 
#' @param B a number of perturbation used for P-value approximation
#' @param rho a numeric vector defining the search grid used in SMMAT-O for SKAT-O (see the SKAT-O paper for details). Default  = c(0, 0.5, 1)
#' @param return_single Logical parameter indicating whether single variant P-value to be returned or not. Default = FALSE.
#' 
#' @return This function returns a dataframe. The row name is the SNP ID, the first column is the prospective score statistics, the second colum is the retrospective score statistics, the third column is the prospective pvalue and the forth column is the restrospective pvalue
#' 
#' @export
#' 
lsrat_test <- function(srat.est,G, weights='beta',impute.method='mean', GRM = NULL, LD = NULL, tests = c("B", "S", "O", "E", "V", "R"), B = 5000,rho = c(0,0.5,1), return_single = FALSE){

# srat.est = m1;
# G=p0_gen$snp.mat;
#  weights='beta';
#  impute.method='mean';
#  GRM = NULL;
#  LD = NULL;
#   tests = c("B", "S", "O", "E", "V", "R")
#  B = 5000;
#  rho = c(0, 0.5,1)

  #value specification
  m<-srat.est$m;time<-srat.est$time;Y = srat.est$Y
  cluster.id<-srat.est$cluster.id;
  snp.names<-colnames(G); family = srat.est$family;
  N = nrow(srat.est$Y);X_1 = srat.est$X; tau = srat.est$tau
  V = srat.est$V; V.inv = srat.est$V.inv; res = srat.est$raw_res;
  mu = srat.est$mu; P1 = srat.est$P1; P2 = srat.est$P2; trans_res = srat.est$trans_res;
  # nullgee = srat.est$nullgee
#symmetrify GRM
	if(!is.null(GRM)){
      GRM = cov2cor(GRM)
    }
#impute G
  G<-as.matrix(G);
  G[G==9]<-NA
  N_MISS<-sum(is.na(G))
  if(N_MISS>0){
    msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G)/ncol(G))
    warning(msg,call.=F)
    G<-snp_impute(G)
  }


var.G<-t(t(G)-colMeans(G));

#center.G projected on X
  colnames(X_1)[1] = "intercept"
  if (is.null(rownames(X_1))){
    rownames(X_1) = paste0('ID', time[,1])
  }
    X_base = do.call(rbind, by(X_1, list(rownames(X_1)), 
                  FUN=function(x) head(x, 1)))

  # calculate center.G based by projecting on X
  # time_cov = grep('Time', colnames(X_1))
  # X_base = as.matrix(X_base[unique(rownames(X_1)), -time_cov])
  # center.G  <- (diag(nrow(X_base))- X_base%*%solve(t(X_base)%*%X_base)%*%t(X_base))%*%G

  center.G = var.G

  SNP.list<-colnames(G)
  MAF<-colMeans(as.matrix(G[,colMeans(var.G^2)!=0]), na.rm=T)/2;
  G<-as.matrix(G[,colMeans(var.G^2)!=0]); center.G = center.G[,colMeans(var.G^2)!=0]

  SNP.name<-SNP.list[colMeans(center.G^2)!=0]
  if(length(weights)>1){Z<-t(t(G)*weights[colMeans(center.G^2)!=0])}
  if(length(weights)==1){
   # toremove.idx = which(MAF>0.05 | MAF<0.001)
   # center.G = center.G[,-toremove.idx]
   # G = G[,-toremove.idx]
   # MAF = MAF[-toremove.idx]   
   weights = rep(0,ncol(G));
   weights<-dbeta(MAF,1,25); 
   Z<-t(t(G)*weights)
   # 09/25/19
   # Z <- t(t(center.G)*weights)
  }
Z<-as.matrix(Z)

#calculate Sigma_G
Sigma_G = cov(center.G)
# grm = cor(t(center.G))

#WG
# weights<-dbeta(maf,1,25);WG<-t(t(G)*weights)
Z0<-as.matrix(Z[match(time[,1],cluster.id),])
p<-ncol(Z0); 
Z<-Z0-P2%*%(P1%*%Z0);

#Prospective V
cluster.id<-unique(time[,1]);m<-length(cluster.id);dimI<-ncol(Z); n.rep<-as.numeric(table(time[,1]));n.total<-1;A = rep(0, m)

    ZPZ<-matrix(0, p, p)

  if(class(srat.est)=="LSRAT_gee"){
    for (i in 1:m)
    {
      ni<-n.rep[i]
      index<-n.total:(n.total+ni-1);
      n.total<-n.total + ni
      Z.i<-Z[index,];if (ni==1) {Z.i<-t(Z.i)}
      Vi.inv<- V.inv[[i]]; res.i = res[index];covy.i<-res.i%*%t(res.i);
      if (family$family =="gaussian"){
        ZPZ<-ZPZ+t(Z.i)%*%Vi.inv%*%covy.i%*%Vi.inv%*%Z.i
      }else{
        mu.i = mu[index]
        if(length(index)>1){delta <- diag((mu.i*(1-mu.i))) }else{delta <- mu.i*(1-mu.i)}  
        ZPZ <- ZPZ+t(Z.i)%*%delta%*%Vi.inv%*%covy.i%*%Vi.inv%*%delta%*%Z.i;
      }
      A[i] = sum(trans_res[index])
    } 
}else{
      for (i in 1:m)
    {
      ni<-n.rep[i]
      index<-n.total:(n.total+ni-1);
      n.total<-n.total + ni
      Z.i<-Z[index,];if (ni==1) {Z.i<-t(Z.i)}
      Vi.inv<- V.inv[[i]]; 
      ZPZ<-ZPZ+t(Z.i)%*%Vi.inv%*%Z.i;
      A[i] = sum(trans_res[index])
    }
}
  


# retrospective V
W = diag(weights)
WGW = W%*%Sigma_G%*%W
if(is.null(GRM)){
  V.retro = c(t(A)%*%A)*WGW
}else{

  V.retro = c(t(A)%*%GRM%*%A)*WGW}


vec.m = rep(1, length(weights))
Q = diag(ncol(W))-WGW%*%vec.m%*%t(vec.m)/sum(WGW)


score.B<-lsrat.perturb(trans_res,Y,time,Z, B, V.inv, family,Q)
score<-score.B$score;B.score<-t(score.B$B.score); score.e  <- score.B$score.e
Q1<-sum(score^2);Q2<-sum(score)^2;lambda_bd = score.B$lambda_bd; Q1.e = sum(score.e^2);
cat("Q1: ", Q1, ",Q2: ", Q2, "\n")
B.Q1<-apply(B.score^2,1,sum);B.Q2<-apply(B.score,1,sum)^2

Burden <- "B" %in% tests
SKAT <- "S" %in% tests
ACAT  <- "V" %in% tests

SKATO <- "O" %in% tests
SMMAT <- "E" %in% tests
ACATR <- "R" %in% tests

#SRAT-B  
if(Burden|SMMAT|SKATO|ACATR){
  p.burden.B = Get.p(as.matrix(Q2), as.matrix(B.Q2));
  p.burden.pro<-Get.p.eigen.burden(Q2, sum(ZPZ/m));
  p.burden.retro <- Get.p.eigen.burden(Q2, sum(V.retro/m));
}
#SRAT-S
if(SKAT|ACATR){
  p.skat.B = Get.p(Q1, B.Q1);
  p.skat.pro  <-  .quad_pval(Q1, ZPZ/m)
  p.skat.retro <- .quad_pval(Q1, V.retro/m)
}
#SRAT-O

if(SKATO){
  # re <- Get.p.eigen.skato(Q1, Q2, ZPZ/m, rho=rho, p.burden.pro)
  re  <- Get.p.eigen.skato2(score, ZPZ/m, rho =rho)
  rho.skato.pro = re$minp.rho
  p.skato.pro = re$p
  # re <- Get.p.eigen.skato(Q1, Q2, V.retro/m, rho=rho,p.burden.retro)
  re <- Get.p.eigen.skato2(score, V.retro/m, rho = rho)
  p.skato.retro = re$p
  rho.skato.retro = re$minp.rho
}


#SRAT-E

if(SMMAT){
  p.smmat.pro = .smmat_pvalue(score, ZPZ/m, p.burden.pro)
  p.smmat.retro = .smmat_pvalue_retro(Q1.e, Q%*%V.retro/m, p.burden.retro)
}


# ACAT Test
if(ACAT|ACATR){
  #create LBRAT est type

  Y.res = srat.est$raw_res;mu = srat.est$mu
  rare.index = which(MAF*m < 10)
  if(class(srat.est)=="LSRAT_gee"){
    #connstruct GEE test
      lbrat.est = list(Y=Y,time=time,X=X_1,cluster.id=cluster.id,m=m, est.type = "GEE", P1=P1, P2=P2, V = V, V.inv = V.inv, Y.res=Y.res, mu = mu, family = family)

    }else{
    # construct GLMM null model
      # lbrat.est = list(Y=Y,time=time,X=X_1,cluster.id=cluster.id,m=m, est.type = "GLMM", P1=P1, P2=P2, V = V, V.inv = V.inv, Y.res=trans_res*tau[3], mu = mu, family = family)
      lbrat.est = list(Y=Y,time=time,X=X_1,cluster.id=cluster.id,m=m, est.type = "GLMM", P1=P1, P2=P2, V = V, V.inv = V.inv, Y.res=trans_res, mu = mu, family = family)
    }

  if(length(rare.index)>0){
        single.P = lbrat_test(lbrat.est, G[,-rare.index]);}else{
        single.P = lbrat_test(lbrat.est, G); 
        }

  single.pro = single.P$pval.pro; single.retro = single.P$pval.retro

  if(length(rare.index>0)){
  burden.pro = Get.p.eigen.burden(sum(score[rare.index])^2, sum(ZPZ[rare.index, rare.index]/m));
  burden.retro = Get.p.eigen.burden(sum(score[rare.index])^2, sum(V.retro[rare.index, rare.index]/m))

  weights.acat = dbeta(MAF[-rare.index],1,25)*(MAF[-rare.index]*(1-MAF[-rare.index]))^0.5
  maf.rare = mean(MAF[rare.index])
  weights.acat = c(dbeta(maf.rare, 1, 25)*sqrt(maf.rare*(1-maf.rare)),weights.acat)
  pvalue.vec.pro = c(burden.pro,single.pro);
  pvalue.vec.retro = c(burden.retro,single.retro);
  }else{
    weights.acat = dbeta(MAF,1,25)*(MAF*(1-MAF))^0.5
    pvalue.vec.pro = single.pro
    pvalue.vec.retro = single.retro
  }



  acat.pro = aca.pvalue(pvalue.vec.pro, weights.acat)
  acat.retro = aca.pvalue(pvalue.vec.retro, weights.acat)
}


if(ACATR){
  acatR.pro <- aca.pvalue(c(p.skat.pro, p.burden.pro, acat.pro), rep(1,3))
  acatR.retro <- aca.pvalue(c(p.skat.retro, p.burden.retro, acat.retro), rep(1,3))
}



# prepare results
results = vector('list', length(tests) )
names(results) = tests
for(i in 1:length(tests)){
#variant-set tests
  if(tests[[i]]=="B") {
    results[[i]] = c(p.burden.B, p.burden.pro, p.burden.retro)
    names(results[[i]]) = c('perturb', 'pro', 'retro')}
  if(tests[[i]]=="S") {
    results[[i]] = c(p.skat.B, p.skat.pro, p.skat.retro)
    names(results[[i]]) = c('perturb', 'pro', 'retro')}

  if(tests[[i]] =="V") {
    if(!return_single){
      results[[i]] = c(acat.pro, acat.retro)
      names(results[[i]]) = c('pro', 'retro')
    }else{
      results[[i]] <- list(acat.pro = acat.pro, acat.retro = acat.retro, single.P = single.P)}}
#omnibus tests
  if(tests[[i]] =="O") {
    results[[i]] = c(rho.skato.pro, p.skato.pro, rho.skato.retro, p.skato.retro)
    names(results[[i]]) = c('rho.pro', 'pro', 'rho.retro', 'retro' )}
  if(tests[[i]] =="E") {
    results[[i]] = c(p.smmat.pro, p.smmat.retro)
    names(results[[i]]) = c('pro', 'retro')}

  if(tests[[i]] == "R"){
    results[[i]] = c(acatR.pro, acatR.retro)
    names(results[[i]]) = c('pro', 'retro')
  }

  }

results$n.marker = length(MAF)

return(results)

}



lsrat.perturb<-function(trans_res,Y,time,Z.I,B,V.inv, family, Q)
{

  Y.res <- trans_res
  N<-nrow(time);
  cluster.id<-unique(time[,1]);m<-length(cluster.id);dimI<-ncol(Z.I)
  ##Generalized score test
  score.array<-score.array.e <- matrix(0,dimI,m)
  n.total<-1;n.rep<-as.numeric(table(time[,1])); lambda_bd = 0;
  for (i in 1:m)
  {
    ni<-n.rep[i]
    index<-n.total:(n.total+ni-1);
    n.total<-n.total + ni
    Z.Ii<-Z.I[index,];if (ni==1) {Z.Ii<-t(Z.Ii)}
    score.array[,i]<-t(Z.Ii)%*%Y.res[index]/sqrt(m)
  }
  score<-as.matrix(apply(score.array,1,sum))
  score.e <- as.matrix(apply(Q%*%score.array,1,sum))
  B.coef<-matrix(rbinom(m*B,1,0.5)*2-1,m,B)
  B.score<-score.array%*%B.coef
  return(list(score=score,B.score=B.score, lambda_bd = lambda_bd, score.e = score.e))#Q: test statistic; M: covariance matrix
}


aca.pvalue  <- function(pvalue, weights){
  rm.idx = which(is.na(pvalue)|pvalue>(1-10^-5));
  if(length(rm.idx)>0){
  pvalue = pvalue[-rm.idx]; weights = weights[-rm.idx]}
  T_aca <- c(weights^2)%*%tan((0.5-pvalue)*pi)
  # cat("T_aca: ", T_aca, '\n')
  p.aca = 0.5-atan(T_aca/sum(weights^2))/pi
}
