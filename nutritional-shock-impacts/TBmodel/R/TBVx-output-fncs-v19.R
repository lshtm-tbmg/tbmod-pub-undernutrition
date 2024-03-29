#output.to.array=function(Y,parms){
#  years = Y[[1]]
#  n     = length(years)
#  idxs  = 1:n
#  dimZ  = c(n,parms$nVXa, parms$nSES, parms$nRISK, parms$nHIV,parms$nTB,parms$nAGES)
#  ndim  = length(dimZ)
#  Z     = array(0,dim=dimZ)
#  tb    = 1:parms$nTB
#  if(ndim==7){
#    for (i in idxs){
#      for (vxa in 1:parms$nVXa)
#        for (ses in 1:parms$nSES)
#          for (risk in 1:parms$nRISK)
#            for (hiv in 1:parms$nHIV){
#              Z[i,vxa,ses,risk,hiv,,]=Y[[2]][[i]][calc.index.by.nr(parms,vxa,ses,risk,hiv,tb),]  
#            }
#    }
#    dimnames(Z)=list(YEAR=years,VXa=names(parms$VXa), SES=names(parms$SES), RISK=names(parms$RISK), HIV=names(parms$HIV), TB=names(parms$TB), AGE=parms$AGES)
#  }else{
#    stop("unsupported number of dimensions in output")
#  }
#  Z
#}

#generate.filename = function(fparams,s,path="."){
#  s1 = unlist(strsplit(fparams$xmlfile,"/",fixed=T))
#  s2 = unlist(strsplit(s1[length(s1)],".xml",fixed=T))
#  s3 = unlist(strsplit(fparams$xmlrunfile,"/",fixed=T))
#  s4 = unlist(strsplit(s3[length(s3)],".xml",fixed=T))
#  s5 = as.character(as.integer(seconds(now())) %% 100000)
#  md5run = digest(fparams$run.params,algo="sha1")
#  md5mdl = digest(fparams,algo="sha1")
#  sig  = paste(substring(md5mdl,nchar(md5mdl)-4),"_",substring(md5run,nchar(md5run)-4),sep="")
#  paste(path,"/output/",s2,"_",s4,"_",sig,"_",s5,s,".txt",sep="")
#}

#' a helper function called by generate.single.dim.flows.output() to generate a data frame with incidence output for a specific dimension
#' 
#' @param t           vector of output time points
#' @param inci.dim    list of matrices of various types of incidence (TBp, Tm etc) split out into dY, dY.in, dY.out for all time points
#' @param dim         the dimension to generate incidence output for
#' @param subjectname subject e.g. "Tm" "TBp" etc
#' @param matname     matrix name i.e. "dY.in" or "dY.out"
#' @param params      initialized model parameters
#' @param cols        set of cols of dimensions with state names by dimension
#' 
#' 
#' @return a data table with incidence output for a specific subject by year and age group 
#' @export
generate.inci.output.df=function(t,inci.dim,dim,subjectname,matname,params,cols){
  assert_that(!is.na(subjectname),msg="subjectname argument in generate.inci.output.df cannot be NA")
  assert_that(matname=="dY.in" | matname=="dY.out", msg = "matname argument in generate.inci.output.df should be dY.in or dY.out")
  
  #if (is.na(subjectname)){
  #  if (matname=="dY.in"){mname="in"}else if(matname=="dY.out"){mname="out"} # else{mname=NA}
  #  flow = rep(mname,nrows)
  #  nrows = nrow(inci.dim[[1]])
  #  matname = NA
  #}else{
  nrows = nrow(inci.dim[[1]][[subjectname]][[matname]])
  #}
  dim   = rep(dim,nrows)
  subject = rep(subjectname,nrows)
  if (matname=="dY.in"){mname="in"}else if(matname=="dY.out"){mname="out"} # else{mname=NA}
  flow = rep(mname,nrows)
  cols = as.data.frame(cbind(cols,dim,subject,flow),stringsAsFactors=T)
  output = params$run.params$output
  sel = which(t %in% output$years)
  assert_that(length(sel)>0,msg="no overlap between years specified for incidence output and detailed output")
  basic_df = data.table(year=rep(0,nrows))
  basic_df = cbind(basic_df,cols)
  df = data.table()
  for (i in sel){
    year = rep(t[i],nrows)
    M    = NULL
    if (sum(is.na(output$age.from))==0){
      if (is.na(subjectname)){
        if(all(output$age.from==as.integer(colnames(inci.dim[[i]])))){
          M=inci.dim[[i]]
        }else{
          M = aggregate.by.age.groups(inci.dim[[i]],output$age.from, sumcols=F, avg = F)
        }
      }else{
        if(length(output$age.from)==length(as.integer(colnames(inci.dim[[i]][[subjectname]][[matname]]))) 
           && all(output$age.from==as.integer(colnames(inci.dim[[i]][[subjectname]][[matname]])))){
          M = inci.dim[[i]][[subjectname]][[matname]]
        }else{
          M = aggregate.by.age.groups(inci.dim[[i]][[subjectname]][[matname]],output$age.from, sumcols=F, avg = F)
        }
      }
    }else{
      if (is.na(subjectname)){
        M = inci.dim[[i]]
      }else{
        M = inci.dim[[i]][[subjectname]][[matname]]
      }
    }  
    colnames(M)=as.integer(colnames(M))
    rownames(M)=rownames(inci.dim[[i]][[1]][[1]])
    rownames_basic_df = paste(as.character(basic_df$VXa),as.character(basic_df$SES),as.character(basic_df$RISK),
                              as.character(basic_df$HIV),as.character(basic_df$TB),sep="-")
    checked_OK = all(rownames(M)==rownames_basic_df)
    assert_that(checked_OK,msg="FATAL error in flow output ; composite dim names don't match inci output")
    if (!checked_OK){
      modlog(level='FATAL', msg="FATAL error in flow output ; composite dim names don't match inci output")
    }
    local_df = cbind(basic_df,M)
    local_df$year = year
    df = rbind(df,local_df)
  }
  y = melt.data.table(df, measure.vars = 10:ncol(df),variable.name = "age_from", value.name = "value", variable.factor = F)
  y$age_from = as.integer(y$age_from)
  
  z = cbind(y[,1:(ncol(y)-1)],age_thru=as.integer(rep(0,nrow(y))),value=y$value)
  agesfrom   = sort(unique(z$age_from))
  agesthru   = c(agesfrom[2:length(agesfrom)],100)-1
  for (i in seq_along(agesfrom)){
    sel = z$age_from==agesfrom[i]
    z[sel,'age_thru']=agesthru[i]
  }    
  z  
}

#' generates prevalence output 
#' 
#' @param t      vector of output time points
#' @param state  list of matrices of state variables x age groups by time points
#' @param params initialized model parameters
#' 
#' 
#' @return a data table with prevalence output by year and age group 
#' @export
generate.prevalence.output=function(t,state,params){
  VXa  = calc.names.for.dim(params,"VXa")
  SES  = calc.names.for.dim(params,"SES")
  RISK = calc.names.for.dim(params,"RISK")
  HIV  = calc.names.for.dim(params,"HIV")
  TB   = calc.names.for.dim(params,"TB")
  cols = cbind(VXa,SES,RISK,HIV,TB)
  z = as.data.table(generate.prev.output.df(t,state,params,cols))
  output = rbind(z,z[,.(age_from=min(age_from),age_thru=max(age_thru),value=sum(value)),by=c(names(z)[1:6])])
  output = cbind(country=factor(params$run.params$countrycode,params$run.params$countrycode),output)
  sel = output$year %% 1 == 0
  output$year[sel] = output$year[sel] - 1e-3
  output
}

#' helper function called by generate.prevalence.output() to generate prevalence output as a data frame 
#' 
#' @param t      vector of output time points
#' @param state  list of matrices of state variables x age groups by time points
#' @param params initialized model parameters
#' @param cols   names of states by individual dimension
#' 
#' 
#' @return a data frame with prevalence output
#' @export
generate.prev.output.df=function(t,state,params,cols){
  nrows = nrow(state[[1]])
  output = params$run.params$output
  sel = which(t %in% output$years)
  assert_that(length(sel)>0,msg="no overlap between years specified for incidence output and detailed output")
  basic_df = data.frame(year=rep(0,nrows))
  basic_df = cbind(basic_df,as.data.frame(cols,stringsAsFactors = T))
  df = data.frame()
  for (i in sel){
    year = rep(t[i],nrows)
    M    = NULL
    if (sum(is.na(output$age.from))==0){
      M = aggregate.by.age.groups(state[[i]],output$age.from, sumcols=F, avg = F)
    }else{
      M = state[[i]]
    }  
    colnames(M)=paste("A",colnames(M),sep="")
    rownames(M)=rownames(state[[i]])
    rownames_basic_df = paste(as.character(basic_df$VXa),as.character(basic_df$SES),as.character(basic_df$RISK),
                              as.character(basic_df$HIV),as.character(basic_df$TB),sep="-")
    checked_OK = all(rownames(M)==rownames_basic_df)
    assert_that(checked_OK,msg="FATAL error in flow output ; composite dim names don't match inci output")
    if (!checked_OK){
      modlog(level='FATAL', msg="FATAL error in flow output ; composite dim names don't match inci output")
    }
    local_df = cbind(basic_df,M)
    local_df$year = year
    df = rbind(df,local_df)
  }
  y = as.data.frame(data.table::melt(setDT(df), measure.vars = patterns("^A\\d+$"),variable.name = "age_from", value.name = "value"))
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
#' generates the flows output for a single dimension 
#' 
#' @param dim.inci incidence output for a single dimension (result from incidence.from.model.run())
#' @param params   initialized model parameters
#' @param dim      dimension name e.g. "TB"
#' 
#' 
#' @return a data frame with flows output
#' @export
generate.single.dim.flows.output=function(t,inci.dim,params,dim=NA){
  VXa  = calc.names.for.dim(params,"VXa")
  SES  = calc.names.for.dim(params,"SES")
  RISK = calc.names.for.dim(params,"RISK")
  HIV  = calc.names.for.dim(params,"HIV")
  TB   = calc.names.for.dim(params,"TB")
  cols = cbind(VXa,SES,RISK,HIV,TB)
  subjects = names(inci.dim[[1]])
  matnames = c("dY.in","dY.out")     
  z    = data.table()
  for (subject in subjects){
    for (matname in matnames){
      if (!is.null(inci.dim[[1]][[subject]])){
        df = generate.inci.output.df(t,inci.dim,dim,subject,matname,params,cols)
        z = rbind(z,df)
      }
    }
  }
  # setDT(z)
  # add 0 thru 99 age group
  minagefrom = min(z$age_from)
  maxagethru = max(z$age_thru)
  z = rbind(z,z[,.(age_from=minagefrom,age_thru=maxagethru,value=sum(value)),by=c(names(z)[1:9])])
  z = cbind(country=factor(params$run.params$countrycode,params$run.params$countrycode),z)
  sel = z$year %% 1 == 0
  z$year[sel] = z$year[sel] - 1e-3
  z
}

#' generates the output for all flows
#' 
#' @param inci   incidence output (result from incidence.from.model.run())
#' @param params initialized model parameters
#' 
#' @return a data frame with flows output
#' @export
generate.all.flows.output=function(inci,params){
  if (!is.null(params$run.params$output)){
    out = list(prev=NULL,TB=NULL,HIV=NULL,VXa=NULL,DEM=NULL)
    for (dim in names(inci)[-(1:2)]){
      if (!is.null(inci[[dim]])){
        gc()
        out[[dim]]=generate.single.dim.flows.output(t=inci$t,inci.dim=inci[[dim]],params,dim=dim)
      }    
    }
    result = out$TB
    if (!is.null(out$HIV)){
      result = rbind(result,out$HIV)
    }
    if (!is.null(out$VXa)){
      result = rbind(result,out$VXa)
    }
    if (!is.null(out$RISK)){
      result = rbind(result,out$RISK)
    }
    if (!is.null(out$DEM)){
      result = rbind(result,out$DEM)
    }
    return(result)
  }
  NULL
}

#' Writes model output to a file
#' 
#' @param data data to be written (data.table or data.frame)
#' @param filename filename
#' @param output.format output format: either 'txt' (default) for tab delimited text, 
#' 'parquet' for parquet and 'fst' for fst
#' 
#' @return the modified list A after adding the result of fn to elements dY, dY.in and dY.out 
#' @export
write.output=function(data,filename, output.format='txt'){
  if (output.format=='parquet'){
    filename = gsub("txt$","parquet",filename)
    write_parquet(data,filename)
  }else if (output.format=='fst'){
    filename = gsub("txt$","fst",filename)
    write.fst(data,path=filename, compress=100)
  }else{
    write.table(data,file=filename,row.names=F,sep="\t",quote=F)
  }  
  print(filename)
}

#' Simple helper function to add elements dY, dY.in and dY.out to an existing list with the same elements 
#' 
#' @param A  a list with elements dY, dY.in and dY.out
#' @param fn a function that produces a list with elements dY, dY.in and dY.out
#' 
#' @return the modified list A after adding the result of fn to elements dY, dY.in and dY.out 
#' @export
add.dY=function(A=NULL,fn=NULL,include=F){
  if (include){
    res      = fn
    if (is.null(A)){
      A = list(dY = 0.*res$dY, dY.in = 0.*res$dY, dY.out = 0.*res$dY)
    }
    A$dY     = A$dY + res$dY
    A$dY.in  = A$dY.in  + res$dY.in
    A$dY.out = A$dY.out + res$dY.out
  }
  A
}

#' Regenerates incidence output from model output 
#' 
#' @param out    model output from run.deSolve() i.e. a list of matrices (state variables by age group) as a function of time
#' @param params environment with initialized model parameters
#' 
#' @return a list with elements t, PREV, and depending on output options e.g. TB, HIV, VXa with 
#' t:    a vector of time points 
#' PREV: a list of matrices by time point with prevalence
#' TB:   a list of incidence elements by time point with TB progression (TBp), TB transmission(Tm) and the various TB treatment matrices (e.g. TBtr_init)
#'       each incidence element contains dY, dY.in and dY.out
#' HIV, VXa, etc: cf. TB     
#' 
#' @export
incidence.from.model.run=function(out,params){
  output = params$run.params$output
  sel    = which(out$t %in% output$years)
  result = list(t=out$t[sel])
  for (i in seq_along(result$t)){
    out.index           = sel[i]
    result$PREV[[i]]    = out$state[[out.index]]
    if (params$contacts$defined){
      params$contacts$M=update.contact.matrix(params, colSums(result$PREV[[i]][params$ALIVE,]))
    }
    df.sel = output$options[output$options$dim=="TB",]
    result$TB[[i]] = incidence.from.Y.matrix(t=out$t[out.index],out$state[[out.index]],params,sel.output=df.sel)
    deaths = result$TB[[i]]$TBdeaths
    result$TB[[i]]$TBdeaths = NULL
    if (params$DIMLENGTHS["HIV"]>1 & any(as.logical(output$options[output$options$dim=="HIV",2:5]))){
      df.sel = output$options[output$options$dim=="HIV",]
      result$HIV[[i]] = incidence.from.Y.matrix(t=out$t[out.index],out$state[[out.index]],params,sel.output=df.sel)
      deaths = deaths + result$HIV[[i]]$HIVdeaths
      result$HIV[[i]]$HIVdeaths = NULL
    }    
    if (params$DIMLENGTHS["VXa"]>1 & any(as.logical(output$options[output$options$dim=="VXa",2:5]))){
      df.sel = output$options[output$options$dim=="VXa",]
      result$VXa[[i]] = incidence.from.Y.matrix(t=out$t[out.index],out$state[[out.index]],params,sel.output=df.sel)
    }
    if (params$DIMLENGTHS["RISK"]>1 & any(as.logical(output$options[output$options$dim=="RISK",2:5]))){
      df.sel = output$options[output$options$dim=="RISK",]
      result$RISK[[i]] = incidence.from.Y.matrix(t=out$t[out.index],out$state[[out.index]],params,sel.output=df.sel)
    }
    
  }
  modlog(level="DEBUG",msg="successful call of incidence.from.model.run()")
  result
}

#' Helper function called by incidence.from.model.run()
#' 
#' @param Y a matrix with state variables
#' @param p environment with initialized model parameters
#' 
#' @return a list with incidence elements such as Tm, TBp, TBdeaths, TBtr_init, etc each with dY, dY.in and dY.out
#' 
#' @export
incidence.from.Y.matrix=function(t,Y,p,sel.output=NA){

  if (sum(is.na(sel.output))>0)
    return(NULL)

  if (sel.output$dim=="TB"){
     out = list()
     out$Tm  = add.dY(A=NULL,fn=derivs.Tm.in.out(t,Y,p),include=sel.output$transmission) 
     out$TBp = add.dY(A=NULL,fn=matmul.by.age.group.in.out(p$TBp,Y),include=sel.output$progression)
     out$TBdeaths = out$Tm$dY.in[p$DEAD,] + out$TBp$dY.in[p$DEAD,]
     for (i in seq_along(p$TBtr)){
       tparams = p$TBtr[[i]]$timed.parameters
       if (!is.null(tparams))
         mult = get(names(tparams),envir=tparams)(t)
       else
         mult = 1
       name = paste0("TBtr_",names(p$TBtr)[i]) 
       out[[name]] = add.dY(A=NULL,fn=matmul.by.age.group.in.out(p$TBtr[[i]],Y,mult),include=sel.output$treatment)
       out$TBdeaths = out$TBdeaths + out[[name]]$dY.in[p$DEAD,]
     }
     out$Xi  = NULL
     for (i in seq_along(p$inci$TB)){
       q   = p$inci$TB[[i]]
       if (q$dim==sel.output$dim){
         out$Xi = add.dY(out$Xi,fn=derivs.inci.in.out(t,Y,p,q,abs(t %% 1)<1e-4),include=sel.output$incidence)
       }
     }
     if (!is.null(out$Xi)){
       out$TBdeaths = out$TBdeaths + out$Xi$dY.in[p$DEAD,]
     }
     return(out)
  }else{
    Tiname = paste0(sel.output$dim,"p")
    out = list()
    out[[Tiname]] = add.dY(A=NULL,fn=matmul.by.age.group.in.out(p[[Tiname]],Y=Y),include=sel.output$progression)
    if (sel.output$dim=="HIV"){
      out$HIVdeaths = out[[Tiname]]$dY.in[p$DEAD,]
      for (i in seq_along(p$HIVtr)){
        name = paste0("HIVtr_",names(p$HIVtr)[i]) 
        tparams = p$HIVtr[[i]]$timed.parameters
        if (!is.null(tparams))
          mult = get(names(tparams),envir=tparams)(t)
        else
          mult = 1
        out[[name]] = add.dY(A=NULL,fn=matmul.by.age.group.in.out(p$HIVtr[[i]],Y, mult),include=sel.output$treatment)
        out$HIVdeaths = out$HIVdeaths + out[[name]]$dY.in[p$DEAD,]
      }
    }
    if (sel.output$dim=="RISK"){
      for (i in seq_along(p$RISKtr)){
        name = paste0("RISKtr_",names(p$RISKtr)[i]) 
        tparams = p$RISKtr[[i]]$timed.parameters
        if (!is.null(tparams))
          mult = get(names(tparams),envir=tparams)(t)
        else
          mult = 1
        out[[name]] = add.dY(A=NULL,fn=matmul.by.age.group.in.out(p$RISKtr[[i]],Y, mult),include=sel.output$treatment)
      }
    }
    Xiname = paste0(sel.output$dim,"Xi")
    out[[Xiname]]  = NULL
    for (i in seq_along(p$inci[[sel.output$dim]])){
      q   = p$inci[[sel.output$dim]][[i]]
      if (q$dim==sel.output$dim){
        out[[Xiname]] = add.dY(A=out[[Xiname]],fn=derivs.inci.in.out(t,Y,p,q,abs(t %% 1)<1e-4),include=sel.output$incidence)
        if (sel.output$dim=="HIV"){
          out$HIVdeaths = out$HIVdeaths + out[[Xiname]]$dY.in[p$DEAD,]
        }
      }
    }
    return(out)
  }    
}


