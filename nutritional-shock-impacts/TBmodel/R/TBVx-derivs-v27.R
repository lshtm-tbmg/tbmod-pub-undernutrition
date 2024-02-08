# Nov 26, 2022. Model version 3.2.8.9

#' Adds vaccination results in case vaccination is once per year ; these data will be part of the output
#'
#' @param t current time 
#' @param p an environment with initialized model parameters
#' @param dvacc an environment with initialized model parameters
#' @param rownamesZ an environment with initialized model parameters
#' @param subject   refers to one of the VXa incidence files e.g. VXaXi_1 etc
#' 
#' @return updated p$vaccinated
#' @export
add.vaccination.results = function(t,p,dvacc,rownamesZ,subject="VXa"){
  ages.from = p$run.params$output$age.from
  M = aggregate.by.age.groups(dvacc,ages.from, sumcols=F, avg = F)
  colnames(M)=paste0("A",colnames(M))
  VXa  = calc.names.for.dim(p,"VXa")
  SES  = calc.names.for.dim(p,"SES")
  RISK = calc.names.for.dim(p,"RISK")
  HIV  = calc.names.for.dim(p,"HIV")
  TB   = calc.names.for.dim(p,"TB")
  assert_that(all(paste(VXa,SES,RISK,HIV,TB,sep="-") == rownamesZ),msg="problem in setting up df for vx reporting: rownames not correct")
  Min   = M ; Min[M<0]=0. ; Mout = M ; Mout[M>0]=0.
  dfin  = cbind(year=t,VXa=VXa,SES=SES,RISK=RISK,HIV=HIV,TB=TB,dim="VXa",subject=subject,flow="in" ,as.data.frame(Min) )
  dfout = cbind(year=t,VXa=VXa,SES=SES,RISK=RISK,HIV=HIV,TB=TB,dim="VXa",subject=subject,flow="out",as.data.frame(Mout))
  result = rbind(dfin,dfout)
  if (is.null(p$vaccinated)){
    return(result)    
  }else{
    return(rbind(p$vaccinated,result))
  }
}
#' Handles aging, adding newborns, adjusting population composition, updating the contacts matrices
#'
#' \code{aging.and.births.and.update.contact.matrix.event()} is called once a year, at the start of a new year
#' 
#' aging.and.births.and.update.contact.matrix.event() takes care of:
#' \describe{
#' \item{\code{aging}:}{moving (part of) the contents of the state matrix \code{Y} (after converting the vector \code{X}) to the next higher age group }
#' \item{\code{births}:}{adding newborns}
#' \item{\code{fiting demography}:}{adjusting the population composition to UNPOP data (baseline) or reusing recorded population adjustments (intervention)}
#' \item{\code{rescaling}:}{rescaling the population (usually in 1950)}
#' \item{\code{updating the contacts matrices}:}{i.e. balancing the contacts matrix to ensure equal number of contacts between ego - alter age groups and vice versa which is necessary if population composition changes}
#' \item{\code{yearly incidence}:}{adding incident cases if incidence is defined as a yearly event rather than a rate}
#' }
#' 
#' @param t current time 
#' @param X a \code{double[nVXaSESRISKHIVTB x nAGES]} vector of state variables
#' @param p an environment with initialized model parameters
#' 
#' @return a \code{double[nVXaSESRISKHIVTB x nAGES]} vector of adjusted state variables
#' @export
aging.and.births.and.update.contact.matrix.event <- function(t,X,p){
  # cat("t=",t,"\n")
  if (t>p$run.params$simulation$years[1]){
    events.at = unique(trunc(p$run.params$simulation$years))[-1]
    if (t==events.at[1] & !p$intervention){
      p$popadj = data.frame(matrix(0,length(events.at),p$nAGES))
      colnames(p$popadj)=p$AGES
      rownames(p$popadj)=events.at
    }
    if (t==events.at[1] & p$intervention & length(p$inci$VXa)>0){
      p$vaccinated = NULL
    }
    Y = Y.matrix.from.X.vector(X,p)

    if (t>=p$run.params$rescale.pop.yr & !p$run.params$rescaled.pop){
      Y = scale.alive.population(Y,p)
      p$run.params$rescaled.pop=T
    }
    # 1. age population (see Schenzle et al.) ; we use deltaY to not confuse with a rate (dY)
    Y[!p$ALIVE,]       = 0. 
    deltaY             = t(tcrossprod(p$aging.matrix, Y))
    # deltaY[!p$ALIVE,]  = 0. # the dead do not age
    deltaY[,1]         = deltaY[,1] + p$birthrate(t)*sum(colSums(Y))*p$at.birth 
    Y                  = Y + deltaY
    # 2. update the population to match UNPOP
    if (p$intervention & t>=p$intervention.start){
      Z = t(p$popadjrateintv[p$popadjfnintv(t),] * t(Y))
    }else{
      unpop    = p$unpopdem[p$unpopdemfn(t),]
      popbyage = colSums(Y)
      if (any(popbyage==0)){
        modlog(level='ERROR',paste0("at yearly event at t = ",t," empty age groups in population:",colnames(Y)[zeros]))
        assert_that(!any(popbyage==0),msg="zeros in popbyage")
      }
      if (!p$intervention){
        p$popadj[as.integer(rownames(p$popadj))==t,] = unpop/popbyage
      }
      Z = t(unpop/popbyage * t(Y)) 
      if (any(is.nan(Z))){
        modlog(level='ERROR',paste0("at yearly event at t = ",t," NaNs in Z"))
      }
    }
    
    deltaZ = 0 * Z
    offset = 1e-4 # should be copied from some other location ....
    for (i in seq_along(p$inci$RISK)){
      if (p$inci$RISK[[i]]$inci.once.per.year){
        drisk = derivs.inci(t+offset,Z,p,p$inci$RISK[[i]],T)
        deltaZ = deltaZ + drisk
      }
    }
    for (i in seq_along(p$inci$HIV)){
      if (p$inci$HIV[[i]]$inci.once.per.year){
        dhiv = derivs.inci(t+offset,Z,p,p$inci$HIV[[i]],T)
        deltaZ = deltaZ + dhiv
      }
    }
    Z = Z + deltaZ
    deltaZ = 0 * deltaZ
    for (i in seq_along(p$inci$TB)){
      if (p$inci$TB[[i]]$inci.once.per.year){
        dtb = derivs.inci(t+offset,Z,p,p$inci$TB[[i]],T)
        deltaZ = deltaZ + dtb
      }
    }
    Z = Z + deltaZ
    deltaZ = 0 * deltaZ
    for (i in seq_along(p$inci$VXa)){
      if (p$inci$VXa[[i]]$inci.once.per.year & t>=p$intervention.start & p$inci$VXa[[i]]$inci.fn(t+offset)>0){
        dvacc  = derivs.inci(t+offset,Z,p,p$inci$VXa[[i]],T)
        colnames(dvacc)=p$AGES
        p$vaccinated = add.vaccination.results(t,p,dvacc,rownames(Z),subject=paste0("VXaXi_",i))
        deltaZ = deltaZ + dvacc
      }
    }
    Z = Z + deltaZ
    if (p$contacts$defined){
      #if (any(colSums(Z[p$ALIVE,])==0)){
      #  cat("t=",t," 0's in Z[p$ALIVE,] in yearly event ...\n")
      #}
      p$contacts$M=update.contact.matrix(p, colSums(Z[p$ALIVE,]))
    }

    if (sum(is.nan(Z))>0){
      modlog(level='FATAL',msg=paste0("at yearly event at t = ",t," NaNs in state variables ; reduce euler / rk2 / rk4 time step"))
    }else{ 
      neg = p$run.params$num.int$min.value.for.state.var
      if (sum(Z<neg)>0){
        idx = which(Z<neg)
        r   = idx %% nrow(Z)
        c   = (idx %/% nrow(Z))+1 
        for (i in seq_along(idx)){
          modlog(level='ERROR',paste0("at yearly event at t = ",t," state variables < ",neg," :",
          " state = ",rownames(Z)[r[i]]," ; age group = ",colnames(Z)[c[i]],"; value = ",signif(Z[r[i],c[i]])," ; please reconsider parameter values, incidence data or setting min.value.for.state.var to a more negative value"))
        }
        if (!(is.null(p$paths$parameters) | is.na(p$paths$parameters))){
          if (!dir.exists(file.path(p$paths$output.dir,"weird"))){
            dir.create(file.path(p$paths$output.dir, "weird"))
          }
          outfname  = paste0(create.filename(paste0(p$paths$output.dir,"/weird"), p, runtype = "weird"),".csv")
          df = data.frame(unique.name = model$modified.params$unique.name, value=unlist(model$modified.params$mean))
          write.csv(df,outfname,row.names = F)
        }
        assert_that(sum(Z<neg)==0,msg = paste("state variables <",neg))
        return(0*X)
      }
    }
    return(as.vector(Z))
  }else{
    return(X)
  }
}

#' Calculate the force of infection for Mtb transmission
#'
#' \code{foi()} calculates the force of infection for Mtb transmission taking into account the age dependent contact matrices
#' 
#' @param t current time 
#' @param Y a \code{double[nVXaSESRISKHIVTB x nAGES]} matrix of state variables
#' @param parms an environment with initialized model parameters
#' 
#' @return a \code{double [nAGES]} of the force of infection by age group
#' @export
foi=function(t,Y,parms){
  tparams = parms$Tm$timed.parameters
  if (!is.null(tparams)){
    mult = get(names(tparams),envir=tparams)(t)
  }else{
    mult = 1
  }
  I_ = matmul.by.age.group(parms$Infc, Y)
  i_ = colSums(I_)/colSums(Y[parms$ALIVE,])
  i_[is.nan(i_) | is.infinite(i_)]=0.
  if (parms$contacts$defined){
    return(parms$DAYSPERYEAR * mult * as.numeric(parms$contacts$M %*% i_)) # return the foi ; foi is a column vector of foi by age 
    # NOTE added 23 July 2023 : code is correct including the code to adjust the contact matrix ; see the update.contact.matrix() function
    # NOTE: this may seem incorrect as originally the contactees are in the columns of the contact matrix read
    #       however adjusting the contact matrix takes care of this ... the new contact matrix M depends on M plus t(M) anyway
    #       and the adjustment takes care of balancing i.e. the number of individuals in a row age group (i.e. susc) multiplied with an
    #       off-diagonal contact rate and the number of individuals in that column balances that same number mirrored wrt the diagonal
    #       effectively, the rows now contain the contactees (susceptibles) and therefore the foi is calculated correctly
  }else{
    avg.i_  = sum(I_)/sum(Y[parms$ALIVE,]) # i.e. the average fraction infectious
    avg.i_[is.nan(avg.i_) | is.infinite(avg.i_)]=0.
    i_      = rep(avg.i_,ncol(Y))
    return(parms$DAYSPERYEAR * mult * as.numeric(i_)) # return the foi ; foi is a column vector with identical average fraction infectious times mult and DAYSPERYEAR i.e. contact rate = 1 / day
  }
}

#' Calculates the derivative dY/dt of matrix Y due to Mtb transmission split into negative and positive derivativs (outflow and inflow)
#'
#' \code{derivs.Tm()} calculates the contribution of Mtb transmission to the derivatives of \code{Y} at time \code{t} with parameters \code{parms} 
#' Within \code{derivs.Tm()} the function \code{foi()} is called to calculate the force of infection.  
#' 
#' @param t current time 
#' @param Y a double[nVXaSESRISKHIVTB x nAGES] matrix of state variables
#' @param parms an environment with initialized model parameters
#' 
#' @return a list of dY (net flow), dY.in (inflow), dY.out (outflow)
#' @export
derivs.Tm.in.out=function(t,Y,parms){
  susc   = matmul.by.age.group.in.out(parms$Tm, Y)
  foi    = foi(t=t,Y=Y,parms=parms)
  dY.in  = t(t(susc$dY.in ) * foi )
  dY.out = t(t(susc$dY.out) * foi )
  return( list(dY=(dY.in+dY.out),dY.in=dY.in,dY.out=dY.out) )
}

#' Calculates the derivative dY/dt of matrix Y due to Mtb transmission
#'
#' \code{derivs.Tm()} calculates the contribution of Mtb transmission to the derivatives of \code{Y} at time \code{t} with parameters \code{parms} 
#' Within \code{derivs.Tm()} the function \code{foi()} is called to calculate the force of infection.  
#' 
#' @param t current time 
#' @param Y a double[nVXaSESRISKHIVTB x nAGES] matrix of state variables
#' @param parms an environment with initialized model parameters
#' 
#' @return a \code{double [nVXaSESRISKHIVTB x nAGES]} of derivatives of state variables
#' @export
derivs.Tm=function(t,Y,parms){
  susc = matmul.by.age.group(parms$Tm, Y)
  foi  = foi(t=t,Y=Y,parms=parms)
  return(t ( t(susc) * foi ) )
}

#' Calculates the derivative dY/dt of matrix Y due to incidence
#'
#' \code{derivs.inci.common()} calculates the contribution of Mtb transmission to the derivatives of \code{Y} at time \code{t} with parameters \code{parms} 
#' 
#' @param t current time 
#' @param Y a \code{double[nVXaSESRISKHIVTB x nAGES] matrix} of state variables
#' @param p an environment with initialized model parameters
#' @param p.inci an environment derived from options in the XML e.g.
#' 
#' \preformatted{
#' <HIV.incidence>
#'   <incidence.data file="data/HIV-incidence.txt" times="1980,2051" values="lambdaH,lambdaH" proportions="false" denominator="susc"/>
#' </HIV.incidence>
#' }
#' 
#' and from the file that defines the incidence (e.g. “data/HIV-incidence.txt”)
#' 
#' The most important data structure is \code{p.inci$inci.matrix} which is based on the incidence data file.
#' 
#' @return a \code{double [nVXaSESRISKHIVTB x nAGES] matrix} of derivatives of state variables
#' @export
derivs.inci.common=function(t=NA,Y=NA,p=NA,p.inci=NA){ 
  dX=NULL
  x_inci = as.numeric(p.inci$inci[p.inci$inci.fn(t),])
  if (sum(x_inci)>0){
    if (!is.null(p.inci$inci.trend.fn)){
      x_inci = x_inci * p.inci$inci.trend.fn(t)
    }
    X = t(Y)
    dX = 0*X
    
    if (p.inci$inci.target.prevalence){ # in this case x_inci really is the target prevalence
        subtract = rowSums(X[,p.inci$tosel])/rowSums(X[,p.inci$allsel])
        subtract[is.nan(subtract)]=0
        dX[,p.inci$susc] = (x_inci - subtract) * X[,p.inci$susc]
    }else{
      if (!p.inci$inci.proportions){
        if (p.inci$inci.denominator=="all"){
          denom   = rowSums(X[,p.inci$susc])
          denom   = sapply(denom, max, 1e-6)
          scaleby = rowSums(X[,p.inci$alive])/denom
          scaleby = sapply(sapply(scaleby, max, 1.),min,10)
        }else{
          scaleby = 1.0
        }
      }
      else {
        if (p.inci$inci.denominator=="all"){
          denom   = rowSums(X[,p.inci$alive])
        }else{
          denom   = rowSums(X[,p.inci$susc])
        }
        denom   = sapply(denom, max, 1)
        scaleby = sum(denom)/denom
        scaleby = sapply(sapply(scaleby, max, 1.),min,50)
      }
      dX[,p.inci$susc]   = x_inci * X[,p.inci$susc] * scaleby # OK - that is simple
    }
  }
  dX
}  

#' Calculates the derivative dY/dt of matrix Y due to incidence 
#'
#' @param t current time 
#' @param Y a double[nVXaSESRISKHIVTB x nAGES] matrix of state variables
#' @param p an environment with initialized model parameters
#' @param p.inci incidence data structure including incidence matrices etc
#' @param once   logical that should match p.inci$once.per.year
#' 
#' @return a double[nVXaSESRISKHIVTB x nAGES] matrix of derivatives of state variables due to incidence
#' @export
derivs.inci=function(t=NA,Y=NA,p=NA,p.inci=NA,once=F){ 
  if (once==p.inci$inci.once.per.year){
    dX = derivs.inci.common(t=t,Y=Y,p=p,p.inci=p.inci)
    if (!(is.null(dX) | all(dX==0))){
      #return(as.matrix(tcrossprod(p.inci$inci.matrix,dX)))
      return(matrix(Matrix::tcrossprod(x=p.inci$inci.matrix,y=dX),ncol(dX),nrow(dX)))
    }
  }
  return(0.*Y)
}

#' Calculates the derivative dY/dt of matrix Y due to incidence split into negative and positive contributions (outflow and inflow)
#'
#' @param t current time 
#' @param Y a double[nVXaSESRISKHIVTB x nAGES] matrix of state variables
#' @param p an environment with initialized model parameters
#' @param p.inci incidence data structure including incidence matrices etc
#' @param once   logical that should match p.inci$once.per.year
#' 
#' @return a list of \code{double[nVXaSESRISKHIVTB x nAGES]} matrices of derivatives of state variables due to incidence i.e. dY (net), dY.in (inflow) and dY.out (outflow)
#' @export
derivs.inci.in.out=function(t=NA,Y=NA,p=NA,p.inci=NA,once=F){ 
  if (once==p.inci$inci.once.per.year){
    dX = derivs.inci.common(t=t,Y=Y,p=p,p.inci=p.inci)
    if (!is.null(dX)){
      P = Q = p.inci$inci.matrix
      P[p.inci$inci.matrix<0]=0.
      Q[p.inci$inci.matrix>0]=0.
      dY.in  = Matrix::tcrossprod(P,dX)
      dY.out = Matrix::tcrossprod(Q,dX) 
      return(list(dY=as.matrix(dY.in+dY.out),dY.in=as.matrix(dY.in),dY.out = as.matrix(dY.out)))     
    }
  }
  return(list(dY=0.*Y,dY.in=0.*Y,dY.out = 0.*Y))     
}

#' Scales the current population composition and size to the initial population composition and size
#' 
#' @param Y a \code{double [p$nVXaSESRISKHIVTB x p$nAGES]} of state variables of the current population
#' @param parms an environment with initialized model parameters
#' 
#' @return a \code{double [p$nVXaSESRISKHIVTB x p$nAGES]} of state variables of the rescaled population
#' @export
scale.alive.population=function(Y,parms){
  Y[!parms$ALIVE,] = 0.
  ini.age.dist     = apply(parms$y.ini,2,sum)
  current.age.dist = apply(Y,2,sum)
  Y = scale(Y,center=F,scale = current.age.dist / ini.age.dist) 
  Y[is.nan(Y) | is.infinite(Y)]=0 
  Y 
}

#' Convert a vector of state variables to a matrix
#' 
#' \code{Y.matrix.from.X.vector()} converts a vector X of state variables to a matrix Y using parameters p
#' 
#' @param X a double[] vector of state variables
#' @param p an environment with initialized model parameters
#' 
#' @return a \code{double [p$nVXaSESRISKHIVTB x p$nAGES]} of state variables
#' @export
Y.matrix.from.X.vector = function(X,p){
  dim(X) = c(p$nVXaSESRISKHIVTB,p$nAGES)
  rownames(X)=p$VXaSESRISKHIVTB
  colnames(X)=names(p$AGES)
  X  
}

#' Calculate dY/dt (a vector of derivates of state variables)
#' 
#' \code{derivs.deSolve()} calculates \code{dY/dt} [a vector of derivatives of state variables] at timepoint \code{t} using parameters \code{p}
#' 
#' Further details:
#' 
#' \describe{
#' \item{\code{Y = Y.matrix.from.X.vector(X,p)}}{convert vector X to matrix Y}
#' \item{\code{dY = matmul.by.age.group(p$TBHIVptr,Y)}}{Calculate dY due to TB and HIV progression and treatment [only for constant matrices] ; \code{TBHIVptr} has at least \code{TBp} and is optionally combined with \code{HIVp} and \code{TBtr} and \code{HIVtr}}
#' \item{\code{dY = dY + derivs.Tm(t,Y,p)}}{TB transmission}
#' \item{\code{if (!is.null(p$TBtr)) \emph{update dY}}}{ with TB treatment (if not already in TBHIVptr)}
#' \item{\code{for (i in seq_along(p$inci$TB))\emph{update dY}}}{ with TB incidence}
#' \item{\code{if (p$nHIV>1)\emph{update dY}}}{ with HIV: progression, treatment, incidence (except when already in HIVTBptr)}
#' \item{\code{if (p$nRISK>1 & !is.null(p$RISKp))\emph{update dY}}}{with RISK progression}
#' \item{\code{if (p$nSES>1  & !is.null(p$SESp))\emph{update dY}}}{with SES progression}
#' \item{\code{if (p$nVXa>1)\emph{update dY}}}{with VXa: progression, incidence}
#' }
#' 
#' @section A special section:
#' 
#' Text for section A
#' 
#' 
#' @param t current time
#' @param X a double[] vector of state variables
#' @param p an environment with initialized model parameters
#' 
#' @export
#' 
#' @return a \code{list[2]} with two \code{double[]} vectors: \code{dY} and \code{dY.HIVTBdeaths}
#' @export
derivs.deSolve=function(t,X,p){
  
  Y = Y.matrix.from.X.vector(X,p) # * p$ALIVE
  #if (any(colSums(Y)==0)){
  #  cat("t=",t," 0's in colSums(Y) ...\n")
  #}
  ntimesteps <<- ntimesteps+1
  dY = matmul.by.age.group(p$TBHIVptr,Y)                # TB and HIV progression and treatment
                                                        # TBHIVptr has at least TBp and is optionally combined with HIVp and TBtr and HIVtr
  dY = dY + derivs.Tm(t,Y,p)                            # TB transmission
  #if (any(is.nan(dY))){
  #  cat("t=",t," NaNs in dY ...\n")
  #}
  if (!is.null(p$TBtr)){                                # TB treatment (if not already in TBHIVptr)
    for (i in seq_along(p$TBtr)){
      if (!p$TBtr[[i]]$aggregated){
        tparams = p$TBtr[[i]]$timed.parameters
        if (!is.null(tparams)){
          mult = get(names(tparams),envir=tparams)(t)
          if (abs(mult)>1e-6){
            dY = dY + matmul.by.age.group(p$TBtr[[i]],Y,mult=mult)
          }
        }else{
          dY = dY + matmul.by.age.group(p$TBtr[[i]],Y)
        }
      }
    }
  }
  for (i in seq_along(p$inci$TB)){                      # TB incidence
    dY = dY + derivs.inci(t,Y,p,p$inci$TB[[i]])
  }
  if (p$nHIV>1){                                        # HIV: progression, treatment, incidence  
    if (!p$HIVp$aggregated){                            # HIV progression only if not yet in HIVTBptr
      dY = dY + matmul.by.age.group(p$HIVp, Y)
    }
    if (!is.null(p$HIVtr)){
      for (i in seq_along(p$HIVtr)){
        if (!p$HIVtr[[i]]$aggregated){                  # HIV progression only if not yet in HIVTBptr
          tparams = p$HIVtr[[i]]$timed.parameters
          if (!is.null(tparams)){
            mult = get(names(tparams),envir=tparams)(t)
            if (abs(mult)>1e-6){
              dY = dY + matmul.by.age.group(p$HIVtr[[i]],Y,mult=mult)
            }
          }else{
            dY = dY + matmul.by.age.group(p$HIVtr[[i]],Y)
          }
        }
      }
    }
    for (i in seq_along(p$inci$HIV)){
      dY = dY + derivs.inci(t,Y,p,p$inci$HIV[[i]])
    }
  }
  if (p$nRISK>1){
    if (!is.null(p$RISKp)){                   # RISK: progression
      dY = dY + matmul.by.age.group(p$RISKp, Y)
    }
    if (!is.null(p$RISKtr)){
      for (i in seq_along(p$RISKtr)){
        if (T){#!p$HIVtr[[i]]$aggregated){          
          tparams = p$RISKtr[[i]]$timed.parameters
          if (!is.null(tparams)){
            mult = get(names(tparams),envir=tparams)(t)
            if (abs(mult)>1e-6){
              dY = dY + matmul.by.age.group(p$RISKtr[[i]],Y,mult=mult)
            }
          }else{
            dY = dY + matmul.by.age.group(p$RISKtr[[i]],Y)
          }
        }
      }
    }
    
    for (i in seq_along(p$inci$RISK)){        # RISK: incidence
      dY = dY + derivs.inci(t,Y,p,p$inci$RISK[[i]])
    }
  }
  if (p$nSES>1  & !is.null(p$SESp)){
    dY = dY + matmul.by.age.group(p$SESp, Y)
  }
  if (p$nVXa>1){                                       # VXa: progression, incidence
    if(!is.null(p$VXap)){
      dY = dY + matmul.by.age.group(p$VXap, Y)
    }
    for (i in seq_along(p$inci$VXa)){
      #dvac = derivs.inci(t,Y,p,p$inci$VXa[[i]])
      #rownames(dvac)=p$VXaSESRISKHIVTB
      #rows = calc.all.indices.for.dim(p,"VXa")==which(p$VXa == "vac")
      dY = dY + derivs.inci(t,Y,p,p$inci$VXa[[i]])
    }
    
  }
  # assert_that(sum(is.nan(dY))==0,msg=paste0("NaNs in derivs.deSolve @ t=",t))
  assert <<- F # disable special type of assert 
  # NOTE: HIV and TB deaths are included in dY (as negative contributions) already due to progression)
  list(dY=as.numeric(dY), dY.HIVTBdeaths=as.numeric(colSums(dY[p$DEAD,])))
}

#' called from run.model() 
#' 
#' within run.deSolve the deSolve function rk() is called with derivs.deSolve() as the argument that provides a reference to the derivs function that produces the derivatives at a specific time point.
#' 
#' @param params an environment with initialized model parameters
#' 
#' @return a \code{list[3]} of \code{double[] times}, \code{list matrix[nVXaSESRISKHIVTB,nAGES] state}, \code{list vector[nAGES] dHIVTBx}
#' @export 
run.deSolve = function(params = NULL){
  assert <<- T
  ntimesteps <<- 0
  started.at = proc.time()
  if (params$contacts$defined){
    params$contacts$M=update.contact.matrix(params,colSums(params$y.ini[params$ALIVE,]))
  }
  params$run.params$rescaled.pop = F
  rawout = NULL
  for (i in 1:length(params$run.params$num.int)){# NOTE there are either 1 or 2 num.int elements
    if (length(params$run.params$num.int)==1){
      output.times = params$run.params$simulation$years
      event.times  = as.numeric(unique(trunc(output.times))[-1])
    }else if (i==1){
      output.times = params$run.params$simulation$years[params$run.params$simulation$years <=   params$run.params$num.int[[2]]$from.year]
      event.times  = as.numeric(unique(trunc(output.times))[-1])
    }else{
      output.times = params$run.params$simulation$years[params$run.params$simulation$years >=    params$run.params$num.int[[2]]$from.year]
      event.times  = as.numeric(unique(trunc(output.times)))
    }
    numint       = params$run.params$num.int[[i]]
    if (i==1){
      yini = as.vector(params$y.ini)
    }else{
      N = params$nVXaSESRISKHIVTB * params$nAGES
      row = nrow(rawout)
      X = rawout[row,2:(N+1)]
      yini = X
    }
    zz           = rk(y        = yini, 
                      times    = output.times, 
                      func     = derivs.deSolve, 
                      parms    = params, 
                      rtol     = numint$rtol, 
                      atol     = numint$atol, 
                      method   = numint$method,
                      maxsteps = numint$maxsteps,
                      hini     = numint$hini,
                      hmin     = numint$hmin,
                      hmax     = numint$hmax,
                      events   = list(func = aging.and.births.and.update.contact.matrix.event, time = event.times))
    rawout       = rbind(rawout,zz)
  }
  if (ntimesteps>(max(params$run.params$simulation$years)-min(params$run.params$simulation$years))*50){
    modlog(level='WARN',msg=paste0(" ntimesteps = ",ntimesteps," in ",timetaken(started.at)," which may be excessive ; please reconsider parameter values"))
  }else{
    modlog(level='INFO',msg=paste0(" ntimesteps = ",ntimesteps," in ",timetaken(started.at)))
  }
  
  t         = rawout[,1]
  out       = vector("list",length(t)) 
  dHIVTBx   = vector("list",length(t)) 
  nr        = params$nVXaSESRISKHIVTB
  nc        = params$nAGES
  N         = nr*nc
  for (i in seq_along(t)){
    tnow = t[i]
    M = matrix(rawout[i,2:(N+1)],nr,nc) ; colnames(M)=params$AGES ; rownames(M)=params$VXaSESRISKHIVTB
    out[[i]] = M
    range = (N+2):(N+nc+1)
    dHIVTBx[[i]]   = as.vector(rawout[i,range])    ; names(dHIVTBx[[i]]) = params$AGES
  }
  sel = diff(t)!=0.
  list(times=t[sel], state=out[sel], dHIVTBx=dHIVTBx[sel])
}

#' Calculates the vaccination flows from p$vaccinated (the yearly update of once per year vaccinations)
#' 
#' @param p an environment with initialized model parameters
#' 
#' @return a data frame of yearly vaccination flows (if option once.per.year=T) to be included in model output e.g. as result$vaccinated
#' @export 
calc.vac.flows=function(p){
  df = p$vaccinated
  y = as.data.frame(melt(setDT(df), measure.vars = patterns("^A\\d+$"),variable.name = "age_from", value.name = "value"))
  z = cbind(y[,1:(ncol(y)-1)],age_thru=as.integer(rep(0,nrow(y))),value=y[,ncol(y)])
  z$age_from = as.integer(substring(z$age_from,2))
  agesfrom   = unique(z$age_from)
  agesthru   = c(agesfrom[2:length(agesfrom)],100)-1
  for (i in seq_along(agesfrom)){
    sel = z$age_from==agesfrom[i]
    z[sel,'age_thru']=agesthru[i]
  }
  z
}

#' called from run() 
#' 
#' within run.model the run.deSolve() function is called to run the model ; in addition the output of run.deSolve() is processed to produces stocks and flows output and in addition output on deaths and counts
#' 
#' @param model.params an environment with initialized model parameters
#' @param output.flows if F no flows output is generated
#' 
#' @return a \code{list} with  \code{stocks} and subsets of stocks i.e. \code{dead.999}, \code{count.999}, \code{alive.500}, \code{alive.999}
#'         and in addition \code{flows} (if \code{output.flows=T}) and if \code{econ.output}="true" in the XML also     
#'         \code{population}, \code{dHIVTBx}, \code{dBGx}, \code{dPOPadj} 
#'         and if vaccination active also \code{vaccinated}
#' @export 
run.model = function(model.params=NULL, output.flows=T){
  out  = run.deSolve(model.params)
  
  modlog(level="DEBUG",msg="successful result of run.deSolve()")
  options = model.params$run.params$output
  for (j in seq_along(out$times)){
    if (length(model.params$run.params$num.int)==1){
      neg = model.params$run.params$num.int[[1]]$min.value.for.state.var
    }else{
      if (out$times[j]<model.params$run.params$num.int[[2]]$from.year)
        neg = model.params$run.params$num.int[[1]]$min.value.for.state.var
      else
        neg = model.params$run.params$num.int[[2]]$min.value.for.state.var
    }
    Y = out$state[[j]]
    tocorrect = Y<0 & Y>neg
    ntocorrect = sum(tocorrect)
    if (ntocorrect>0){
      modlog(level='WARN',paste(" at t =",out$times[j],ntocorrect,"state variables < 0 and > ",neg," reset to 0 after simulation run and before processing output"))
    }
    Y[Y<0 & Y>neg]=0
    out$state[[j]]=Y
    if (sum(Y<0)>0){
      idx = which(Y<0)
      r   = idx %% nrow(Y)
      c   = (idx %/% nrow(Y))+1 
      for (i in seq_along(idx)){
        modlog(level='ERROR',paste0(" at t = ",out$times[j]," state variables < min.value.for.state.var :",
                                    " state = ",rownames(Y)[r[i]],
                                    " ; age group = ",colnames(Y)[c[i]],
                                    "; value = ",signif(Y[r[i],c[i]])))
      }
    }
  }
  result = list()
  result$stocks = generate.prevalence.output(t=out$times,state=out$state,params=model.params)
  dead = rep(F,model.params$nVXaSESRISKHIVTB)
  for (name in model.params$DIMNAMES){
    dead = dead | grepl("dead$",result$stocks[[name]])
  }
  count = rep(F,model.params$nVXaSESRISKHIVTB)
  for (name in model.params$DIMNAMES){
    count = count | grepl("count$",result$stocks[[name]])
  }
  # split stocks into alive.5, alive.999, dead.999, count.999
  result$dead.999  = result$stocks[ dead & signif(result$stocks$year %% 1, 3) == 0.999,]
  result$count.999 = result$stocks[count & signif(result$stocks$year %% 1, 3) == 0.999,]
  result$alive.500 = result$stocks[!(count|dead) & signif(result$stocks$year %% 1, 1) == 0.5,]
  result$alive.999 = result$stocks[!(count|dead) & signif(result$stocks$year %% 1, 3) == 0.999,]
  gc()
  # if (options$suppress_zeros_stocks) { stocks = stocks[value>1e-6,] }
  modlog(level="DEBUG",msg="successful generation of stocks output from generate.prevalence.output()")
  if (output.flows){
    result$flows  = generate.all.flows.output(incidence.from.model.run(out,model.params),params=model.params)
    gc()
  # if (options$suppress_zeros_flows){ flows = flows[abs(value)>1e-6,] }
    modlog(level="DEBUG",msg="successful generation of flow output from generate.flow.output()")
  }
  if (options$econ.output){ 
    demography.output = new.demography.from.model.run(out,params=model.params) 
    result$population = demography.output$pop
    result$dHIVTBx    = demography.output$dHIVTBx
    result$dBGx       = demography.output$dBGx
    result$dPOPadj    = demography.output$dPOPadj
    #if (model.params$intervention & !is.null(model.params$vaccinated)){
    if (!is.null(model.params$vaccinated)){
      result$vaccinated  = calc.vac.flows(model.params)
    } 
  }
  gc()
  result
}

#' Creates demography output from a model run
#' 
#' @param out    
#' @param params an environment with initialized model parameters
#' 
#' @return a \code{list} with \code{pop} i.e. population composition, \code{dHIVTBx} i.e. number of HIV and TB deaths per year, 
#'         \code{dBGx} number of background deaths per year, \code{dPOPadj} fractional population adjustment rates (fraction of the population per year)
#' 
#' @export 
new.demography.from.model.run=function(out,params){
  sel         = (out$times %% 1) == 0 | (out$times %% 1) == 0.5
  t           = out$times[sel]
  indexes     = (1 : length(out$state))[sel]
  assert_that(sum(sel)>0, msg="no time points at ####.5 or ####.0 found ....")
  Y           = matrix(0,nrow=length(t),ncol=ncol(out$state[[1]]))
  colnames(Y) = as.integer(colnames(out$state[[indexes[1]]]))
  selendyear  = (t %% 1) == 0
  tmod        = t
  tmod[selendyear]=tmod[selendyear]-1e-3
  rownames(Y) = tmod
  country     = params$run.params$countrycode
  dHIVTBx     = Y # HIV and TB deaths (numbers / yr)
  dBGx        = Y # BG  deaths (numbers / yr)
  dPOPadj     = Y # fractional adj rates (fraction of the pop / yr)
  
  if (!params$intervention){  
    for (i in seq_along(tmod)){ 
      Y[i,]         = colSums(out$state[[indexes[i]]][params$ALIVE,])
      dHIVTBx[i,]   = out$dHIVTBx[[indexes[i]]]
      dBGx[i,]      = get.deathrate.allcauses(params,tmod[i]) * Y[i,] - dHIVTBx[i,]
      n = sum(dBGx[i,]<0)
      if (n>0){
        modlog(level='DEBUG',msg=paste0("at t=",t[i]," ",n," age groups with HIV TB deaths exceeding all cause mortality (WHO deathrates) ; capping background deaths to ensure background deaths >= 0"))
      }
      dBGx[i,] = pmax(0,dBGx[i,])
      # if (t[i]>3020){cat("BASELINE: t=",t[i],"background deaths=",sum(dBGx[i,]),"\t","BASELINE: t=",t[i],"HIV TB     deaths=",sum(dHIVTBx[i,]),"\n")}
      if (i>1){
        dPOPadj[i,] = as.numeric(params$popadj[as.integer(rownames(params$popadj))==trunc(tmod[i]),])
      }else{
        dPOPadj[i,] = rep(1,ncol(Y))
      }
      # if (t[i]>2020){ cat("t=",t[i]," pop adj = ",dPOPadj[i,1:8],"\n") }
    }
  }else{
    for (i in seq_along(tmod)){
      Y[i,]         = colSums(out$state[[indexes[i]]][params$ALIVE,])
      dHIVTBx[i,]   = out$dHIVTBx[[indexes[i]]]
      dBGx[i,]      = get.bgdeathrate(params,tmod[i]) * Y[i,]
      # if (t[i]>3020){ cat("INTERVENTION: t=",t[i],"background deaths=",sum(dBGx[i,]),"\t","INTERVENTION: t=",t[i],"HIV TB     deaths=",sum(dHIVTBx[i,]),"\n")}
    }
  }
  pop          = cbind(year=tmod,country=country,as.data.frame(Y,row.names = F))
  dHIVTBx.df   = cbind(year=tmod,country=country,as.data.frame(dHIVTBx,row.names = F))
  dBGx.df      = cbind(year=tmod,country=country,as.data.frame(dBGx,row.names = F))
  if (!params$intervention){
    dPOPadj.df = cbind(year=as.numeric(rownames(params$popadj)),country=country,as.data.frame(params$popadj))
  }else{
    dPOPadj.df = NULL
  }
  modlog(level="DEBUG",msg="successful call of demography.from.model.run()")
  list(pop=pop, dHIVTBx=dHIVTBx.df, dBGx=dBGx.df, dPOPadj=dPOPadj.df)
}