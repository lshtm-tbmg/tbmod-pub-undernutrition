# Nov 27, 2022. Model version 3.2.8.9

#' adjust age groups of incidence file to age groups defined in params environment
#' 
#' initializes VXa, SES, RISK, HIV, or TB progression from the XML specification
#' 
#' @param  p   params environment
#' @param  M   matrix with lower limit of age group in column names
#' @param  compress.1st logical to indicate if incidence in age groups < 5 should be compressed in 1st i.e. death at birth
#' @param  proportions  logical to indicate if numbers per row are proportions (and should add up to 1.0)
#' @param  aggregate    logical to indicate if single year age groups should be aggregated to p$AGES
#' @return a matrix with age groups according to p$AGES
#' @export
adjust.age.groups = function(p=NULL,M=NULL,compress.1st=F, proportions=F, aggregate=T){
  if (proportions){
    assert_that(sum(abs(apply(M,1,sum)-1.)>1e-3)==0,msg="incidence proportions do not add up to 1.0 (before redistribution over age groups)")
  }
  x = matrix(0,nrow=nrow(M),ncol=100)
  xages = 0:99
  colnames(x)=xages
  rownames(x)=rownames(M)
  ages = as.integer(colnames(M))
  nages = length(ages)
  diffs = ages[2:nages]-ages[1:(nages-1)]
  start = 1
  if (compress.1st){
    if (!proportions){
      x[,1] = M[,1]*diffs[1]
    }else{
      x[,1] = M[,1]
    }
    xrange = xages>ages[1] & xages<ages[2]  
    x[,xrange]=0
    start = 2
  }
  for (i in start:nages){
      xrange = xages>=ages[i] 
      if (i < nages)
        xrange = xrange & xages<ages[i+1]  
      x[,xrange] = M[,i]
      if (proportions){
        x[,xrange] = M[,i]/sum(xrange)
      }
  }
  if (proportions){
    assert_that(sum(abs(apply(x,1,sum)-1.)>1e-3)==0,msg="incidence proportions do not add up to 1.0 (after redistribution over age groups)")
  }
  if (aggregate){
    y = aggregate.by.age.groups(x,p$AGES, sumcols = F, avg = !proportions)
    if (proportions){
      assert_that(sum(abs(apply(y,1,sum)-1.)>1e-3)==0,msg="incidence proportions do not add up to 1.0 (after aggregating new age groups)")
    }
    rownames(y) = rownames(M)  
    return(y)
  }else{
    return(x)    
  }
}

#' Updates the contacts matrix to ensure that the number of contacts between A1 and A2 is equal to the number of contacts between A2 and A1
#' 
#' 
#' This function is called once per year at the start of the year in function aging.and.births.and.update.contact.matrix.event()
#' 
#' initializes VXa, SES, RISK, HIV, or TB progression from the XML specification
#' 
#' @param  parms   parameters environment
#' @param  alive.pop current alive population size by age group
#' @return updated contacts matrix
#' @export
update.contact.matrix=function(parms, alive.pop){
  n = length(alive.pop)
  A = matrix(alive.pop,n,n)
  # Checked. The number of contacts between A1 and A2 should be equal to the number of contact between A2 and A1
  # The contact matrix is adjusted to ensure that. The average number of contacts should remain roughly the same for not too large diffs in group size
  M = as.matrix(t(parms$contacts$MplustM  *  alive.pop) / (A+t(A)))
  # NOTE: thought the t(parms$contacts$MplustM  *  alive.pop) may seem incorrect it is definitely correct 
  # in the foi function the foi calculation includes M %*% i_ (with i_ the fraction infectious)
  # and later the multiplication with susc (i.e. numbers resulting from matrix multiplication of transition matrix and pop (Y))
  # means that all works correctly.
  # check with  
  M[is.nan(M)]=0
  assert_that(isSymmetric(M * alive.pop),msg="adjusted contacts matrix multiplied with pop should be symmetric")
  #assert_that(sum(is.nan(M))<1,msg="too many NaN in adjusted contacts matrix")
  M
}

#' Creates a matrix to do aging once per year (see Schenzle et al.) ; called from function aging.and.births.and.update.contact.matrix.event()
#' 
#' 
#' @param  p   parameters environment
#' @return aging matrix
#' @export
create.aging.matrix = function(p){
  M = diag(x=-c(1/(p$AGES[2:p$nAGES]-p$AGES[(2:p$nAGES)-1]),0.))
  for (i in 1:(p$nAGES-1))
    M[i+1,i]=-M[i,i]
  M
} 

#' Creates a list of vectors of indices of age groups (taking into account age groups size)
#' 
#' 
#' @param  p   parameters environment
#' @param  age.groups lower limits of age groups e.g. a named vector with names c(A0, A15) and values c(0, 15)
#' 
#' @return aging matrix
#' @export
get.age.indices=function(p,age.groups){
  age.from.index = age.groups*0
  age.thru.index = age.groups*0
  age.ranges = list()
  for (ag in seq_along(age.groups)){
    if (ag  == length(age.groups)){
      select = p$AGES >= age.groups[ag] & p$AGES < 200
      age.from.index[ag] = min((1:p$nAGES)[select])
      age.thru.index[ag] = max((1:p$nAGES)[select])
    }
    else{
      select = p$AGES >= age.groups[ag] & p$AGES < age.groups[ag+1]
      age.from.index[ag] = min((1:p$nAGES)[select])
      age.thru.index[ag] = max((1:p$nAGES)[select])
    }
    age.ranges[[ag]] = age.from.index[ag]:age.thru.index[ag]
  }
  age.ranges
}

#' Multiples the age specific matrices in subject with the appropriate subset of Y to produce one matrix of the same dimensions as Y 
#' and splits the result in net, positive and negative contributions dependent on the sign of the multiplying age group matrix 
#' 
#' @param  subject an environment from which the parameter matrices are derived
#' @param  Y a matrix with states in rows and ages in columns
#' @param  mult a multiplier to be applied to the resulting matrix
#' 
#' @return a list of dY (net), dY.in (positive contributions) and dY.out (negative contributions)
#' @export
matmul.by.age.group.in.out=function(subject,Y,mult=1){
  X     = subject$parameter.matrices
  nag   = length(subject$age.groups)
  dY    = dY.in = dY.out = 0. * Y
  if (nag==1){
    Z = as.matrix(X[[1]])
    P = Q = Z
    P[Z<0] = 0.
    dY.in       =  P%*%Y
    Q[Z>0] = 0.
    dY.out      =  Q%*%Y
  }else{
    for (i in 1:nag){
      range = subject$age.ranges[[i]]
      Z = as.matrix(X[[i]])
      P = Q = Z
      P[Z<0]    = 0.
      dY.in[,range]  =  P%*%Y[,range]
      Q[Z>0]    = 0.
      dY.out[,range] =  Q%*%Y[,range]
    }
  }
  return(list(dY=mult*(dY.in+dY.out),dY.in=mult*dY.in,dY.out=mult*dY.out))
}

#' Multiples the age specific matrices in subject with the appropriate subset of Y to produce one matrix of the same dimensions as Y 
#' 
#' 
#' @param  subject an environment from which the parameter matrices are derived
#' @param  Y a matrix with states in rows and ages in columns
#' @param  mult a multiplier to be applied to the resulting matrix
#' 
#' @return the result of the multiplication of subsets of Y with the appropriate matrix for an age group
#' @export
matmul.by.age.group=function(subject,Y,mult=1){
  X     = subject$parameter.matrices
  nag   = length(subject$age.groups)
  if (nag==1){
    return(mult*matrix(X[[1]] %*% Y,nrow(Y),ncol(Y)))
  }else{
    dY    = 0*Y
    for (i in 1:nag){
      range = subject$age.ranges[[i]]
      dY[,range] = matrix(X[[i]] %*% Y[,range],nrow(Y),length(range)) 
    }
    return(mult*dY)
  }
}

#' Aggregates columns with single year age groups into age groups defined by the lower age limits
#' with options to 
#' 
#' @param Y a matrix with single year age groups 0..99 in columns
#' @param lower.age.limits new age group lower limits 
#' @param sumcols logical to indicate if the summed columns will be returned
#' @param avg logical to indicate if the average of the rows in an age group will be returned
#' @param weights an optional matrix with weights (of the same dimensions as Y)
#' 
#' @return the result of the aggregation of columns of matrix Y into the appropriate age groups, optionally summed into a vector and optionally averaged and/or weighted
#' @export
aggregate.by.age.groups=function(Y,lower.age.limits, sumcols=NA, avg = NA, weights=NA){
  nag = length(lower.age.limits)
  Y = as.data.frame(Y)
  if (assert){
    # print("assertions are active in sum.by.age.group")
    if (nag<2 || prod(lower.age.limits>=0)!=1 || prod((lower.age.limits[2:nag]-lower.age.limits[1:(nag-1)])>0)!=1)
      stop("nr of age.limits should be > 1, all age limits should be > 0 and strictly increasing")
    if (sum(as.integer(colnames(Y)) ==  0:99)<100)
      stop("expecting a matrix with 0:99 as column names i.e. age groups")
  }
  ages.0.99 = as.integer(colnames(Y))
  x = matrix(0,nrow(Y),nag)
  colnames(x)=lower.age.limits
  for (i in 1:nag){
    range = ages.0.99 >=lower.age.limits[i] 
    if (i<nag)
      range = range & ages.0.99 <lower.age.limits[i+1]
    if (sum(range)>1){
      if (avg){
        if (sum(is.na(weights)>0)){
          x[,i]=rowSums(Y[,range])/sum(range) # range should be a vector of T/F values i.e. the sum equals the number of columns selected
        }else{
          x[,i]=rowSums(Y[,range]*weights[,range])/rowSums(weights[,range])
        }
      }else{
        x[,i]=rowSums(Y[,range]) 
      }
    }else if (sum(range)==1){
      x[,i]=Y[,range]
    }
  }
  if (sumcols)
    return(colSums(x))
  else
    return(x)
}

#' Matrix multiplication of X and Y
#' 
#' @param  X a square parameters matrix for one age group
#' @param  Y a matrix (of state by ages)
#' 
#' @return the matrix result of X and Y
#' @export
matmul.for.one.age.group=function(X,Y){
  return(as.matrix(X %*% Y))
}

#' Reformats the contacts matrix with 5 year age groups A0..A75 (as read from Polymod data) into the appropriate format
#' 
#' @param  parms the initialized parameters environment
#' @param  contactM contacts matrix with 5 year age groups
#' @param  scale logical that determines scaling
#' 
#' @return the contacts matrix in the appropriate format i.e. parms$nAGES x parms$nAGES
#' @export
expand.contact.matrix = function(parms, contactM, scale=T){ # performance is not important as this is done once per parameters initialization
  M = matrix(0,parms$nAGES,parms$nAGES)
  rownames(M)=colnames(M)=names(parms$AGES) # i.e. 0,1,....,70,80
  ages = as.integer(substring(colnames(contactM),2)) # e.g. 0, 5, 10, ....
  # check contacts matrix
  assert_that(length(ages)==16 & sum((ages[2:16]-ages[1:15])==5)==15)
  # and model matrix
  one.yr   =  parms$nAGES == 82 && sum(parms$AGES != c(0:80,90))==0
  five.yr  =  parms$nAGES == 16 && sum(parms$AGES != 5*(0:15))==0
  first.yr.single = parms$nAGES == 17 && sum(parms$AGES != c(0,1,5*(1:15)))==0
  
  assert_that(one.yr || five.yr || first.yr.single, msg="choose either standard model age groups c(0:80,90) or 5*(0:15) or c(0,1,5*(1:15))")
  if (one.yr){
    for (i in seq_along(ages)){
      if (i<length(ages)){
        range1  = parms$AGES >= ages[i] & parms$AGES < ages[i+1]
        range1a = NULL
      }else{
        range1  = parms$AGES >= ages[i] & parms$AGES < ages[i]+5 # i.e. 80  
        range1a = parms$AGES >= 80 
      }
      for (j in seq_along(ages)){
        if (j<length(ages)){
          range2 = parms$AGES >= ages[j] & parms$AGES < ages[j+1]
          range2a= NULL
        }else{
          range2 = parms$AGES >= ages[j] & parms$AGES < ages[j]+5 # i.e. 80  
          range2a= parms$AGES >= 80 
        }
        M[range1,range2]=contactM[i,j] / 5.
        if (!is.null(range1a)){
          M[range1a,range2]=contactM[i,j]
        }
        if (!is.null(range2a)){
          M[range1,range2a]=contactM[i,j]
        }
        if (!is.null(range1a) & !is.null(range2a)){
          M[range1a,range2a]=contactM[i,j]
        }
      }
    }
  }else if (five.yr){
    M=contactM
  }else if (first.yr.single){
    M[3:17,3:17]=contactM[2:16,2:16]
    M[1: 2,3:17]=contactM[1   ,2:16]
    M[2:17,1   ]=contactM[1:16,1]*0.2
    M[2:17,2   ]=contactM[1:16,1]*0.8
    M[1   ,1:2 ]=M[2,1:2]
  }
  M  
}