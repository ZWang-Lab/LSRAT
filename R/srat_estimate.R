#' GEE NULL model estimation
#'
#' This function estimate the parameters and residuals for the NULL model for GEE/LSRAT tests
#'
#' @param y.long Long-formatted phenotype vector 
#' @param time Time covarites matched with phenotype vector
#' @param y.cov Covariate matrix denoting the covariate variables measured at each time
#' @param timecov Logical variable, indicating whether the time fixed effect is estimated 
#' @param corstr String, correlation structure for GEE model, optional values are: 'ar1', 'ind', 'mixture' 
#' @return This function returns a list object with model parameters and residuals of the NULL GEE model 
#' @export
#' 
lsrat_est <- function(y.long, time, y.cov, timecov = TRUE, corstr = "ar1"){

	# determine binary or continous
    if(length(table(y.long))>2){
      family = gaussian(); scale.fix = F
      cat('Phenotype is continuous, fitting identity link......\n')
    }else{
      family = binomial(); scale.fix = T
      cat('Phenotype is dichotomous, fitting logistic link ......\n')

    }

    # determine design matrix include time effect
    Y<-as.matrix(y.long);
    if(timecov == TRUE){
      X=as.matrix(cbind(y.cov, time[,2]-1));
      colnames(X)[ncol(X)]="Time";
    }else 
    {
      X <- as.matrix(y.cov);
    }

    # add intercept
    N<-length(y.long)
    X_1 = cbind(rep(1, nrow(X)),X)
    cluster.id<-unique(time[,1]);m<-length(cluster.id)

	  nullgee0<-geeglm(Y~X, family=family, id=time[,1], corstr ="ar1", scale.fix = scale.fix)

  #transformed Residual
    n.total<-1;n.rep<-as.numeric(table(time[,1]));res<-rep(0, n.total);
    rho = as.numeric(summary(nullgee0)$corr[1]);disper = as.numeric(summary(nullgee0)$dispersion[1]);P1<-matrix(0, ncol(X_1), N); V<-list();V.inv<-list(); mu  <- nullgee0$fitted.values; Y_mu <- Y-mu;
   for(i in 1:m)
    {
      ni <- n.rep[i];index<-n.total:(n.total+ni-1);n.total<-n.total + ni
      Y_mu.i <- Y_mu[index];covy.i<-Y_mu.i%*%t(Y_mu.i); mu.i <- mu[index]
      sigma <- rho^abs(outer(time[index,2], time[index,2],"-"));
      if(nullgee0$family$family=="gaussian"){
        #Gaussian
        Vi <- sigma*disper; Vi.inv<-solve(Vi); 
        if(length(index)>1){
            P1[,index]<-t(X_1[index,])%*%Vi.inv;}else{
            P1[,index]<-t(X_1[index,])/disper;
        } 
        res[index] <- Vi.inv%*%Y_mu.i;
      }else
        #Binomial
      { 
        if(length(index)>1){delta <- diag((mu.i*(1-mu.i))) }else{delta <- mu.i*(1-mu.i)}          
        Vi <- sqrt(delta)%*%sigma%*%sqrt(delta);Vi.inv<-solve(Vi); 
        if(length(index)>1){

            P1[,index]<-t(X_1[index,])%*%delta%*%Vi.inv%*%covy.i%*%Vi.inv%*%delta;}else{
            P1[,index]<-t(X_1[index,])*c(covy.i);
        } 
        res[index] <- delta%*%Vi.inv%*%Y_mu.i;
      }
      
      #save results

      V[[i]]<-Vi;V.inv[[i]]<-Vi.inv;
    }   
  P2<-X_1%*%solve(P1%*%X_1)
	lsrat_null <- list(Y=Y,time=time,X=X_1,cluster.id=cluster.id,m=m,corstr=corstr, P1 = P1, P2 = P2, trans_res = res, raw_res = Y_mu, V = V, V.inv = V.inv, mu = mu,family = family)
  class(lsrat_null) = "LSRAT_gee"
  return(lsrat_null)

}


#' GLMM NULL model estimation
#'
#' This function estimate the parameters and residuals for the NULL model in SMMAT/RSMMAT test
#'
#' @param y.long Long-formatted phenotype vector 
#' @param time Time covarites matched with phenotype vector
#' @param y.cov Covariate matrix denoting the covariate variables measured at each time
#' @param timecov Logical variable, indicating whether the time fixed effect is estimated
#' 
#' @return This function returns a list object with model parameters and residuals of the NULL GLMM model 
#' @export


rsmmat_est <- function(y.long, time, y.cov, timecov = TRUE){

  Y<-as.matrix(y.long);
  if(timecov == TRUE){
    X=as.matrix(cbind(y.cov, time[,2]-1));
    colnames(X)[ncol(X)]="Time";
  }else 
  {
    X <- as.matrix(y.cov);
  }


  if(length(table(y.long))>2){
    family = gaussian(); scale.fix = F
    cat('Phenotype is continuous, fitting identity link......\n')
  }else{
    family = binomial(); scale.fix = T
    cat('Phenotype is dichotomous, fitting logistic link ......\n')
  }

  N<-length(y.long);X_1 = cbind(rep(1, nrow(X)),X)
  cluster.id<-unique(time[,1]);m<-length(cluster.id);
  subj = time[,1]; rep = time[,2];

  if(family$family=="gaussian"){
    # nullglmm0<-lmer(Y ~ X+(rep|subj));
    nullglmm0 <- lmer(Y ~ X + (1|subj) + (1|rep))

    # extract parameters
    var.df = as.data.frame(VarCorr(nullglmm0))
    tau<-c(var.df$vcov[2],var.df$vcov[1],var.df$vcov[3]); var_e = var.df$vcov[3];

  }else{
    nullglmm0 <- glmer(Y ~ X + (1|subj) + (1|rep), family = family)
    # try1 = try( glmer(Y ~ X+(1|subj)+(1|rep), family = family))
    # if (class(try1) != "try-error"){nullglmm0 <-try1}else{nullglmm0 = glmer(Y ~ X + (1|subj/rep), family = family, nAGQ = 0)}
    var.df = as.data.frame(VarCorr(nullglmm0))
    tau<-c(var.df$vcov[2],var.df$vcov[1])
  }

  beta<-nullglmm0@beta;mu<-nullglmm0@resp$mu;Y_mu <- (Y-mu);rho = 0.7; tau = c(tau, rho)
  # obtain V, V.inv, P1, P2
  n.total<-1;n.rep<-as.numeric(table(time[,1]))
  V<-list();V.inv<-list();P1<-matrix(0, ncol(X_1), N); res<-rep(0, n.total);
    for (i in 1:m)
      {
        ni<-n.rep[i];index<-n.total:(n.total+ni-1);n.total<-n.total + ni
        mu.i<-mu[index];
        v1<-rho^abs(outer(time[index,2], time[index,2],"-"));
        v2<-matrix(1,length(index), length(index));
        if(family$family=="gaussian"){
          #gaussian
          v3 <- diag(lengths(mu.i))
          res[index] <- Y_mu[index]/var_e;
          if(ni >1){ # deal with subj with 1 obs
            sigma<-tau[1]*v1+tau[2]*v2+var_e*v3;}else{
            sigma<- tau[1] + tau[2] + var_e
            }
        
        }else{
          #binomial
        
        res[index] <- Y_mu[index]

        res[index] <- Y_mu[index]
  
        if(ni >1){ # deal with subj with 1 obs
            v3<-diag((mu.i*(1-mu.i))^-1)
            sigma<-tau[1]*v1+tau[2]*v2+v3;}else{
            v3 <- (mu.i*(1-mu.i))^-1
            sigma<- tau[1] + tau[2] + c(v3)
            }
        }


        Vi<-sigma; Vi.inv<-solve(Vi);  
        V[[i]]<-Vi;V.inv[[i]]<-Vi.inv;

        # calculate P1 and P2

        if(ni>1){
          P1[,index]<-t(X_1[index,])%*%Vi.inv;
        }else{
          P1[,index] <- t(X_1[index,])*c(Vi.inv)
        }
      }


  P2<-X_1%*%solve(P1%*%X_1)
  rsmmat_est  <-  list(Y=Y,time=time,X=X_1,cluster.id=cluster.id,m=m,coef = beta, tau = tau, est.type = "GLMM", P1=P1, P2=P2, V = V, V.inv = V.inv, mu = mu, nullglmm = nullglmm0, trans_res = res, raw_res = Y_mu, family = family)
  class(rsmmat_est) = "RSMMAT_glmm"
  return(rsmmat_est)
}


