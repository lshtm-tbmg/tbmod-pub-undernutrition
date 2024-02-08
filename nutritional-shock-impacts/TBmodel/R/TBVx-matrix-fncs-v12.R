#' A helper function called from create.incidence.matrix()
#' 
#' @param M a square matrix
#' @param by a number to roll down (by>0) or up (by<0) a copy of the matrix 'from' -1 values and (by default) negate those values
#' @param negate default T so the 'to' values will be +1's
#' 
#' @return a double matrix of the same dimensions as M with +1s added in the appropriate positions
#' 
roll.matrix.by.rows = function(M,by,negate=T){
  
  # standard: 1:n
  # (by+1):n
  # 1:(n-by)
  assert_that(by!=0,msg="by==0 in roll.matrix.by.rows ???")
  n = nrow(M)
  P = 0*M
  if (by>0){ # roll down
    assert_that(sum(abs(M[(n-by+1):n,]))==0,msg="ERROR in roll.matrix.by.rows in roll.down!!!")
    P[(by+1):n ,] = M[1:(n-by),] 
  }else if (by<0){ # roll up
    by = -by
    assert_that(sum(abs(M[1 :by,]))==0,msg="ERROR in roll.matrix.by.rows in roll up!!!")
    P[1 :(n-by),] = M[(by+1):n ,]  
  }
  if (negate){
    return(M-P)
  }else{
    return(M+P)
  }
}

colSums = function(M){
  if (is.null(dim(M)))
    return(sum(M))
  else
    return(base::colSums(M))
}


#' Creates an incidence matrix (with -1 to indicate 'from' state and +1 to indicate 'to' states)
#' 
#' @param p an environment with initialized model parameters
#' @param dim a string with the name of the dimension e.g. "VXa" or "HIV"
#' @param from the name of the from state e.g. "HIV-"
#' @param to the name(s) of the to state(s) e.g. "HIVu1" and "HIVcount"
#' 
#' @return a double p$nVXaSESRISKHIVTB x p$nVXaSESRISKHIVTB matrix to be used to implement state transitions due to (VXa, HIV) incidence just with -1 and +1 in the appropriate positions
#' 
create.incidence.matrix=function(p,dim,from,to,sel=NULL){ # this function should work also for VXa (or any other dim)
  # code could be improved ....
  assert_that(length(from)==1,msg="specify one 'from' state in incidence file")
  assert_that(length(to)<=2,  msg="specify at most two 'to' states in incidence file of which the 2nd ends in 'count'")
  
  index      =  which(p$DIMNAMES==dim)
  #Q          =  matrix(0,p$nVXaSESRISKHIVTB,p$nVXaSESRISKHIVTB)
  from.idx   =  p$ALIVE & grepl(from,p$inci.dim.names[[index]]) 
  to.idx     =  p$ALIVE & grepl(to[1],  p$inci.dim.names[[index]])
  if (!is.null(sel)){
    from.idx = from.idx & sel
    to.idx   = to.idx & sel
  }
  x = numeric(p$nVXaSESRISKHIVTB)
  x[from.idx]=-1
  M = diag(x)
  d = which(to.idx)[1]-which(from.idx)[1]
  M = roll.matrix.by.rows(M,by=d)
  assert_that(  sum(colSums(M[,from.idx])!=0)==0 & 
                sum(colSums(abs(M[,from.idx]))!=2)==0,
                msg=paste("Error in creating incidence submatrix for ",dim," from=",from," to=",to,sep=""))
  if (length(to)>1){
    assert_that(grepl("count$",to[2]),msg="2nd to should be a state ending in 'count'")
    to.idx     =  p$COUNT & grepl(to[2],  p$inci.dim.names[[index]])
    if (!is.null(sel)){
      from.idx = from.idx & sel
      to.idx   = to.idx & sel
    }
    x = numeric(p$nVXaSESRISKHIVTB)
    x[from.idx]=-1
    Q = diag(x)
    d = which(to.idx)[1]-which(from.idx)[1]
    Q = roll.matrix.by.rows(Q,by=d)
    diag(Q)=0
    assert_that(sum(colSums(Q[,from.idx])!=1)==0,msg=paste("Error in creating incidence count submatrix for ",dim," from=",from," to=",to,sep=""))
    M = M + Q
    assert_that(sum(colSums(abs(M[,from.idx]))!=3)==0,msg=paste("Error in creating composite incidence matrix for ",dim," from=",from," to=",to,sep=""))
  }
  Matrix(M)
}

# Pulls out the bare names of the stages of a specific dimension 
# 
# @param p an environment with initialized model parameters
# @param dimname a string with the name of the dimension e.g. "VXa" or "HIV"
# 
# @return a vector of length p$nVXaSESRISKHIVTB with names of the stages of the specified dimension
# 
# dim.names.for.all=function(p,dimname){
#  dimindex = which(p$DIMNAMES==dimname)
#  names    = rep(NA,p$nVXaSESRISKHIVTB)
#  indices  = calc.all.indices.for.dim(p,dimname)
#  for (i in 1:p$DIMLENGTHS[dimindex])
#    names[indices==i]=p$DIMNAMESLIST[[dimindex]][i]
#  names
#}

#' Calculates the indices (e.g. 1..10 for HIV) of a specific dimension for all stages of the dimension 
#' 
#' @param p an environment with initialized model parameters
#' @param dimname a string with the name of the dimension e.g. "VXa" or "HIV"
#' 
#' @return a vector of length p$nVXaSESRISKHIVTB with indices within the specified dimension
#' 
calc.all.indices.for.dim=function(p,dimname){
  ndim = length(p$DIMLENGTHS)
  v = (1:p$nVXaSESRISKHIVTB)-1
  switch(dimname,
         "VXa" =  return(1+((v %/% prod(p$DIMLENGTHS[2:ndim])) %% p$nVXa )),
         "SES" =  return(1+((v %/% prod(p$DIMLENGTHS[3:ndim])) %% p$nSES )),
         "RISK"=  return(1+((v %/% prod(p$DIMLENGTHS[4:ndim])) %% p$nRISK)),
         "HIV" =  return(1+((v %/% prod(p$DIMLENGTHS[5:ndim])) %% p$nHIV )),
         "TB"  =  return(1+( v %% p$nTB ))
  )
  NULL
}

#' Pulls out the bare names of the stages of a specific dimension from p$VXaSESRISKHIVTB
#' 
#' @param p an environment with initialized model parameters
#' @param dimname a string with the name of the dimension e.g. "VXa" or "HIV"
#' 
#' @return a vector of length p$nVXaSESRISKHIVTB with names of the stages of the specified dimension
#' 
calc.names.for.dim=function(p,dimname){
  names    = p$DIMNAMESLIST[[which(p$DIMNAMES==dimname)]]
  v        = calc.all.indices.for.dim(p,dimname)
  names[v]
}  

#' Creates the indices for one full dimension and a subset of single other dimensions
#' 
#' Called from create.full.parameter.matrix.for.age.group()
#' 
#' @param p an environment with initialized model parameters
#' @param a named list of dimensions all of length 1 except one of full length
#' 
#' @return the indices of a subset of p$VXaSESRISKHIVTB to be used to create a small transition / progression matrix 
calc.small.matrix.indices.for.one.full.dim.and.specific.single.dim.indices.from.list=function(p,list){
  (list$VXa-1) *p$nSES*p$nRISK*p$nHIV*p$nTB + 
    (list$SES-1) *p$nRISK*p$nHIV*p$nTB +
    (list$RISK-1)*p$nHIV*p$nTB + 
    (list$HIV-1) *p$nTB + list$TB
}

#' 
#' 
#' Called from create.full.parameter.matrix.for.age.group()
#' 
#' @param M square matrix (e.g. transition / progression / treatment matrix) 
#' @param from the name of the 'from' state
#' @param to   the name of the 'to' state
#' @param rate the (additional) transition from M[from,from] to M[to,from] i.e. to be subtracted from M[from,from] and to be added to M[to,from]
#'  
#' @return the modified matrix M 
transition = function(M=NULL, from=NA, to=NA, rate=NULL){
  from=as.character(from)
  to  =as.character(to)
  if (!grepl("count$",to)){
    M[from,from]=M[from,from]-rate
  }
  M[  to,from]=M[  to,from]+rate
  M
}

