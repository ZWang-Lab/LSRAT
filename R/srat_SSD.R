#########################
#   Testing PLINK files using SSD file 
#########################
lsrat.SSD.OneSet_SetIndex = function(SSD.INFO, SetIndex, lsrat_est.obj){

  id1<-which(SSD.INFO$SetInfo$SetIndex == SetIndex)
  if(length(id1) == 0){
    MSG<-sprintf("Error: cannot find set index [%d] from SSD!", SetIndex)
    stop(MSG)
  }
  SetID<-SSD.INFO$SetInfo$SetID[id1]

  try1<-try(Get_Genotypes_SSD(SSD.INFO, SetIndex, is_ID=T),silent = TRUE)
  if(class(try1) != "try-error"){
    G<-try1
    Is.Error<-FALSE
  } else {
    err.msg<-geterrmessage()
    msg<-sprintf("Error to get genotypes of %s: %s",SetID, err.msg)
    stop(msg)
  }
  re<-lsrat_test(lsrat_est.obj, G)

  return(re)
}

#'LSRAT test or RSMMAT test using SSD format files
#'
#' 
#'
#' @param SSD.INFO SSD format information file, output of function "Open_SSD". The genome wide scan are run set by set.
#' @param lsrat_est.obj ouput from lsrat_est() or rsmmat_est()
#' @param ... Other options of the LSRAT or RSMMAT test. Defined as same as in function "lsrat_test()" or "rsmmat_est()"
#' 
#' 
#' 
#' @return reults of the LSRAT test or RSMMAT test. First column contains batchID, second column contains SNP ID, third column concains prospective P-value and forth column contains retrospective P-value
#' 
#' @export 
lsrat.SSD.All = function(SSD.INFO, lsrat_est.obj, ...){
  N.Set<-SSD.INFO$nSets
  OUT.pvalue<-matrix(NA, N.Set, 12)
  OUT.Marker = rep(NA, N.Set)
  Is.Error = TRUE
  for(i in 1:N.Set){
    if(i%%100==0){print(paste0(i," sets finished"))}

    try1 = try(lsrat.SSD.OneSet_SetIndex(SSD.INFO=SSD.INFO, SetIndex=i, lsrat_est.obj=lsrat_est.obj))

    if(class(try1) != "try-error"){
      re<-try1;
      Is.Error<-FALSE
    } else {

      err.msg<-geterrmessage()
      msg<-sprintf("Error to run LSRAT for %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
      warning(msg,call.=FALSE)

    }
    print(paste(i,"-th set",SSD.INFO$SetInfo[i,2],",",SSD.INFO$SetInfo[i,3],"SNPs"))

    if(!Is.Error){
      OUT.pvalue[i, ] <- c(re$B[2:3], re$S[2:3], re$V, re$O[c(2,4)], re$E, re$R)
      print(OUT.pvalue[i,])
      OUT.Marker[i]<-re$n.marker
    }
  }
  out.tbl<-data.frame(SetID=SSD.INFO$SetInfo$SetID, P.value=OUT.pvalue, N.Marker=OUT.Marker)
  if(class(lsrat_est.obj)=="RSMMAT_glmm"){
      colnames(out.tbl)[2:13]<-paste0(rep(c('GLMM.', 'RSMMAT.'), 6), rep(c("B", "S", "C", "O", 'E', "A"), each = 2))}else{
      colnames(out.tbl)[2:13]<-paste0(rep(c('GEE.', 'LSRAT.'), 6), rep(c("B", "S", "C", "O", 'E', "A"), each = 2))     
  }
  out.tbl$n.marker = OUT.Marker
  return(out.tbl)
}





