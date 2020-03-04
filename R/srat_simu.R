#' Simulation for LSRAT and RSMMAT test
#'
#' This function use pre-defined parameters to make the simulation data with longitudinal continuous traits for the LSRAT test (including type I and power test)
#'
#' @param n.sample Numeric, sample size, number of individuals
#' @param n.time Numeric, number of measurements for each individual
#' @param par List, the parameters for the phenotype traits, including covaraites and individual specific time dependent random effects
#' @param time_cov Logical variable, indicating whether time effect is included in phenotypic traits
#' @param snp.count Numeric, number of SNPs in each variant set
#' @param intercept Logical variable, indicating whether intercept is used in phenotypic traits
#' @param power Logical variable, indicating whether the phenotype generated under the null model (type I error) or under the alternative model (power test)
#' 
#' @return A list object is returned to be used as object for association test
#' @export

lsrat_simu.conti<-function(n.sample=5000, n.time=7, par=list(),
    time_cov = TRUE, snp.count = 100, intercept=TRUE, power = FALSE, seed = NULL)
{

    if(missing(par) || length(par)==0 )
    {
        par <- list(b0 = 1.0, b1 = 0.5, btime = 2, 
            sig.a = 0.8, sig.b = 0.8, sig.e = 0.8, rho=0.7,
            coef.causal = .05, positive.ratio = .5, causal.prop = .5);    
    }


        if(par$causal.prop ==0.5){
            par$coef.causal = 0.04
        }
        if(par$causal.prop ==0.2){
            par$coef.causal = 0.06
        }
        if(par$causal.prop ==0.05){
            par$coef.causal = 0.12
    }

    # print(par)
    snp.mat <- simu_snp(n.sample, snp.count, seed = seed);

    if(!power){
 
            phe<- simu.phe.null(n.sample, n.time, par, intercept, time_cov);
    }else{
    
            phe<- simu.phe.power(n.sample, n.time, par, intercept, time_cov, snp.mat)
            snp.mat = phe$snp.mat;
    }

    colnames(phe$y) <- paste("Y", 1:(NCOL(phe$y)), sep="") ;
    rownames(phe$y) <- paste("ID", 1:NROW(phe$y), sep="");

    colnames(phe$cov) <- paste("X", 1:(NCOL(phe$cov)), sep="") ;
    rownames(phe$cov) <- rep(paste("ID", 1:n.sample, sep=""), each = n.time);

    y.long <- c(t(phe$y));
    y.time <- cbind(rep(1:nrow(phe$y), each =n.time), rep(1:n.time, nrow(phe$y)))

    return(list(phe.wide = phe$y, phe.long = y.long, phe.time = y.time, phe.cov.long=phe$cov, snp.mat = snp.mat, h2 = phe$h2))
   
}

#' Simulation for LSRAT and RSMMAT test
#'
#' This function use pre-defined parameters to make the simulation data with longitudinal binary traits for the LSRAT test (including type I and power test)
#'
#' @param n.sample Numeric, sample size, number of individuals
#' @param n.time Numeric, number of measurements for each individual
#' @param par List, the parameters for the phenotype traits, including covaraites and individual specific time dependent random effects
#' @param time_cov Logical variable, indicating whether time effect is included in phenotypic traits
#' @param snp.count Numeric, number of SNPs in each variant set
#' @param intercept Logical variable, indicating whether intercept is used in phenotypic traits
#' #' @param power Logical variable, indicating whether the phenotype generated under the null model (type I error) or under the alternative model (power test)
#' 
#' @return A list object is returned to be used as object for association test
#' @export

lsrat_simu.bi <- function(n.sample=5000, n.time=7, par=list(),
    time_cov = TRUE, snp.count = 100, intercept=TRUE, power = FALSE){

       if(missing(par) || length(par)==0 )
    {
        # par <- list(b0 = -1.0, b1 = 0.5, btime = 0.2, 
        par <- list(b0 = -2.5, b1 = 0.5, btime = 0.2,
            sig.a = 0.8, sig.b = 0.8, sig.e = 0, rho=0.7,
            coef.causal = .1, positive.ratio = 1, causal.prop = 1);    
    } 

        # power parameter for binary trait
        if(par$causal.prop ==0.5){
            par$coef.causal = 0.08
        }
        if(par$causal.prop ==0.2){
            par$coef.causal = 0.12
        }
        if(par$causal.prop ==0.05){
            par$coef.causal = 0.24
        }


    phe<- simu.phe.bi.ascertain(n.sample, n.time, par, intercept, time_cov, snp.count, power)


   return(list(phe.wide = phe$phe.wide, phe.long = phe$phe.long, phe.time = phe$phe.time, phe.cov.long=phe$phe.cov, snp.mat = phe$snp.mat, h2 = phe$h2))

}

simu.bi <- function(n.sample, n.time, par, intercept, time_cov, snp.count, power, seed=NULL){
    if (!power){
    mu  <-  simu.phe.null(n.sample, n.time, par, intercept, time_cov)
    snp.mat <- simu_snp(n.sample, snp.count, seed)
    }else{
    snp.mat <- simu_snp(n.sample, snp.count, seed);
    mu  <-  simu.phe.power(n.sample, n.time, par, intercept, time_cov, snp.mat)
    }

    mu.long  <- c(t(mu$y))
    y.long <- rbinom(n.sample*n.time, 1, inv.logit(mu.long))
    y.wide <- matrix(y.long, nrow = n.sample, ncol = n.time, byrow=T)
    return(list(y=y.wide, cov=mu$cov, r = mu$r, h2 = 0, snp.mat = snp.mat));
}



simu.phe.bi.ascertain <- function(n.sample, n.time, par, intercept, time_cov, snp.count, power){
    sampled.y =c(); sampled.cov = c(); sampled.random.y =c();sampled.snp = c()
    total.samp = 0;

    batch.size = 1; seed = round(sample(10000,1))
    while(total.samp<n.sample){
        phe.temp = simu.bi(n.sample, n.time, par, intercept, time_cov, snp.count, power, seed) 
        y = phe.temp$y; cov = phe.temp$cov;snp.mat = phe.temp$snp.mat

        index.g1 <-which(y[,1] ==1);
        index.g0 <- which(y[,1]==0);
        n.case <- min(length(index.g1), round((n.sample-total.samp)/2))
        index.select.1 <- index.g1[sample(n.case)];
        index.select.0 <- index.g0[sample(n.case)];
        ind.sel<-c(index.select.1, index.select.0)

        sampled.y = rbind(sampled.y, y[ind.sel,])
        sampled.cov = rbind(sampled.cov, cov[covSelector(ind.sel, n.time),])
        sampled.snp = rbind(sampled.snp,snp.mat[ind.sel,])
        total.samp = nrow(sampled.y)
        cat("Perform Acertain Sampling..", min(total.samp/n.sample,1)*100, '%\n')
    }

    cat('finished!\n')
    sampled.y = sampled.y[1:n.sample,]
    sampled.cov = sampled.cov[covSelector(1:n.sample, n.time),]
    sampled.snp <- sampled.snp[1:n.sample,]

    # prepare for output
    colnames(sampled.y) <- paste("Y", 1:n.time, sep="") ;
    rownames(sampled.y) <- paste("ID", 1:n.sample, sep="");
    colnames(sampled.cov) <- paste("X", 1:ncol(sampled.cov), sep="") ;
    rownames(sampled.cov) <- rep(paste("ID", 1:n.sample, sep=""), each = n.time);

    y.long <- c(t(sampled.y));
    y.time <- cbind(rep(1:n.sample, each =n.time), rep(1:n.time, n.sample))


    return(list(phe.wide = sampled.y, phe.long = y.long, phe.time = y.time, phe.cov=sampled.cov, snp.mat = sampled.snp))

}


 # generate random effect
f.simu<-function( n.sample, n.time, par)
{
    ncol <- n.time;
    AR1 <- array(0,dim=c(ncol,ncol));
    for(i in 1:ncol)
    for(j in 1:ncol)
        AR1[i,j] <- par$rho^abs(i-j);

    sigma.b <- par$sig.b^2*AR1;

     
    r <- rnorm(n.sample, 0, par$sig.a ) %*% array(1, dim=c(1,ncol))+
         rmvnorm(n.sample,  rep(0, ncol), sigma.b ) 

    if (par$sig.e!=0){
        sigma.e <- diag(par$sig.e^2, ncol);
        r <- r + rmvnorm(n.sample, rep(0, ncol), sigma.e);
    }
    return(r);
}


simu_snp<-function(n.sample, snp.count, seed)
{
    file.snp.hap1 <- system.file("extdata", "skat-test-1.hap.gz", package="LBRAT");
    file.snp.pos2 <- system.file("extdata", "skat-test-1.pos", package="LSKAT");

    snp.pos <- read.table(file.snp.pos2, header=T);
    snp.hap <- read.table(file.snp.hap1, header=F);

    snp.maxpos <- max(snp.pos$CHROM_POS);
    snp.minpos <- min(snp.pos$CHROM_POS);

 

    if(!is.null(seed)){set.seed(seed)} #only setseed the selection of column (SNPs)
    snp.start <- as.integer(runif(1, min=snp.minpos, max=snp.maxpos-50*1000));
    snp.ends <- snp.start + 50*1000;
    p.sel <- which( snp.pos$CHROM_POS>snp.start & snp.pos$CHROM_POS<snp.ends);
    # message("snp selection id: ", head(p.sel))

    set.seed(NULL) #remove setsedd before select row (subj)
    snp.mat1<-snp.hap[ sample(nrow(snp.hap),n.sample), p.sel];
    snp.mat2 <- snp.hap[ sample(nrow(snp.hap), n.sample), p.sel];
    snp.mat <- abs(snp.mat1-2) + abs(snp.mat2 -2);
    #check var(g) !=0       
    maf <- colMeans(snp.mat)/2;
    m.same <- which( maf==1 |  maf==0 );

    if (length(m.same)>0) snp.mat <- snp.mat[, -m.same, drop=F ];
    snp.mat <-snp.mat[,1:snp.count]
    maf = maf[1:snp.count];
    rare.count = sum(maf<0.01)

    cat('Simulated ', ncol(snp.mat), 'SNPs, ', rare.count, " rare SNPs\n")
    rownames(snp.mat)<-paste("ID", 1:nrow(snp.mat), sep = "");
    colnames(snp.mat)<-paste("SNP", 1:ncol(snp.mat), sep = "")
    return(snp.mat)
}



simu_snp.ld<- function(n.sample, snp.count = 1000, ld = 0.5){
    ld.matrix = (1-ld)*diag(snp.count) + ld*rep(1, snp.count)%*%t(rep(1, snp.count))
    latent_cont1 = rmvnorm(n.sample, rep(0, snp.count), ld.matrix)
    latent_cont2 = rmvnorm(n.sample, rep(0, snp.count), ld.matrix)

    maf.rare = runif(n = floor(snp.count*0.5), 0.01, 0.1 )
    maf.common = runif(n= snp.count-floor(snp.count*0.5), 0.1, 0.5)
    maf = c(maf.rare, maf.common)

    snp.mat1 <-  snp.mat2  <- matrix(0, n.sample, snp.count)
    for(j in 1:snp.count){snp.mat1[,j] = ifelse(pnorm(abs(latent_cont1[,j]), lower.tail = F)<maf[j]/2,1, 0)}
    for(j in 1:snp.count){snp.mat2[,j] = ifelse(pnorm(abs(latent_cont2[,j]), lower.tail = F)<maf[j]/2,1, 0)}
    snp.mat <- snp.mat1 + snp.mat2 ;
    return(snp.mat)
}


simu.phe.null<-function( n.sample, n.time, par, intercept, time_cov){


    # cov.mat <- cbind(rnorm(n.sample*n.time, 0, 1 ));
    cov.mat  <- cbind( rnorm(n.sample*n.time, 0, 1 ), rep(ifelse(runif(n.sample)>0.5, 0, 1), each = n.time));
    y.random = f.simu(n.sample, n.time, par)
    if(intercept){
        y <- y.random + matrix(cbind(1, cov.mat )%*%c( par$b0, par$b1,par$b1), n.sample, n.time, byrow = TRUE)
    }else{
        y <- y.random + cov.mat %*% c(par$b1, par$b1);
    }

    if(time_cov == TRUE){
        time.effect <- rep(1, n.sample)%*%(t(seq(0, n.time-1)*par$btime));
        y = time.effect+y;
    } 

    rownames(cov.mat) <- rep(paste("id", 1:n.sample, sep=""), each = n.time);
    return(list(y=y, cov=cov.mat, r = y.random, h2 = 0));
}



simu.phe.power  <- function(n.sample, n.time, par, intercept, time_cov, snp.mat){
    y.null = simu.phe.null(n.sample, n.time, par, intercept, time_cov)
    y0 = y.null$y; cov.mat = y.null$cov; r = y.null$r; 
    maf <- colMeans(snp.mat)/2;
    maf <- ifelse(maf>=.5, 1-maf, maf)
    nc = ncol(snp.mat);ld <- cor(snp.mat); 
    high_ld_idx = unique(which(ld>=0.5& ld!=1, arr.ind = TRUE)[,1]);
    n.causal = round(par$causal.prop*nc)
    cat(length(high_ld_idx), "SNPs are in high LD (>.5)\n")
    snp.mat.c = as.matrix(scale(snp.mat, scale = FALSE))
    sign.causal <- rep(0, nc);
    #proportion of causal SNP
    if(length(high_ld_idx)>n.causal){
        sign.causal[sample(high_ld_idx, n.causal)] = 1;
    }else{
        sign.causal[c(high_ld_idx, sample(setdiff(1:nc, high_ld_idx), n.causal-length(high_ld_idx)))] = 1;
    }

    #proportion of positive signal
    sign.causal = sign.causal*ifelse(runif(nc)<par$positive.ratio,1, -1)
    snp.effect  <- rep(0, nc)
    snp.index = which(maf>0)
    snp.effect[snp.index]  <- sign.causal[snp.index]*par$coef.causal*abs(log10(maf))
    y.snp  <- apply(snp.mat.c, 1, function(x) sum(x*snp.effect))
    y = y0 + y.snp%*%array(1, dim=c(1,n.time))
    h2 = var(y.snp)/(var(y.snp) + var(y0[,1]))
     cat( "h2=", round(h2,3)*100,"%\n");
    return(list(y = y, cov= cov.mat, h2 = h2, snp.mat = snp.mat))
}


