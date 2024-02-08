# Nov 28, 2022. Model version 3.2.8.9									  

#' Attempts to combine age dependent parameter matrices to a single parameter matrix 
#' 
#' @param  Ti an environment with a list of parameter matrices by age grouping
#' 
#' @return updated environment with updated list of parameter matrices
optimize.params=function(Ti){ # for example the transition matrix ...
  M.list = Ti$parameter.matrices
  n = length(M.list)
  if (n==1)
    return(Ti)
  M = M.list[[1]]
  n.elements = nrow(M)*ncol(M)
  uniq = 0*(1:n)+T
  uniq[1]=F
  for (i in 2:length(M.list)){
    if (is.na(sum(sum(M.list[[i]]==M)==n.elements))){
      print("NA in function optimize.params() ???")
    }            
    if (sum(M.list[[i]]==M)==n.elements)
      uniq[i]=F
    else
      return(Ti)
  }
  if (sum(uniq)==0){
    Ti$parameter.matrices = list(M)
    Ti$age.groups = Ti$age.groups[1]
    Ti$age.from.index = Ti$age.from.index[1]
    Ti$age.thru.index = Ti$age.thru.index[n]
    return(Ti)
  }
}

#' Sets the initial fractions of the population (by state and age group) 
#' 
#' Sets the fractions of people in each state and age group using the XML fraction.at.birth attribute of 
#' the VXa, SES, RISK and HIV dimensions and the TB <seeded.infections> XML element.
#' This distribution is used to initialize the population with a prespecified distribution over the various states by age
#' 
#' @param  p   params environment
#' @param  df   data frame from XML seeded infections for TB
#' 
#' @return a matrix of fractions by state and age group
fractions.by.stage.at.start=function(p, df=NULL){

  assert_that(!is.null(df),msg="error in fraction.by.age.at.start(): data frame with initial TB distribution required ")
  
  rownames = p$VXaSESRISKHIVTB
  colnames = p$DIMNAMES[!p$DIMNAMES=="TB"] # exclude TB
  M = matrix(0,length(rownames),length(colnames),dimnames=list(rownames,colnames))

  for (dim in colnames){
    ind = calc.all.indices.for.dim(p,dim)
    for (i in 1:p$DIMLENGTHS[dim])
      M[,dim][ind==i]=p$ATBIRTH[[dim]][i]
  }

  x = apply(M,1,prod) # we now have a vector of fractions of all dims except TB

  P = matrix(0, nrow = p$nVXaSESRISKHIVTB, ncol = p$nAGES, dimnames = list(p$VXaSESRISKHIVTB,names(p$AGES)))
  for (i in 1:nrow(df)){
    rows = calc.all.indices.for.dim(p,"TB")==which(p$TB == df$stage[i])
    range = p$AGES >= df$age.from[i] & p$AGES <= df$age.thru[i]
    P[rows,range] = df$fraction[i]
  }
  y = x * P
  assert_that(sum(colSums(y)-1.0)<1e-9,msg="something rotten in setting up the at.start matrix .....")
  y 
}

#' Sets the fractions of newborns (by state) 
#' 
#' Sets the fractions of newborns in each state using the XML fraction.at.birth attribute of 
#' the VXa, SES, RISK, HIV and TB dimensions
#' This distribution is used to generate an appropiate distribution of newborns over states
#' 
#' @param  p  params environment
#' 
#' @return a vector of fractions by state
fractions.by.stage.at.birth=function(p){
  rownames = p$VXaSESRISKHIVTB ; colnames = p$DIMNAMES
  M = matrix(0,length(rownames),length(colnames),dimnames=list(rownames,colnames))
  for (dim in p$DIMNAMES){
    ind = calc.all.indices.for.dim(p,dim)
    for (i in 1:p$DIMLENGTHS[dim])
      M[,dim][ind==i]=p$ATBIRTH[[dim]][i]
  }
  x = apply(M,1,prod) # we now have a vector of fractions of all dims 
  x
}

replace.NAs.by=function(nodes,attrib,replace=0,first=NA,n=NA){
  v = as.numeric(xml_attr(nodes,attrib))
  if (sum(is.na(v))==n){
    v = rep(0.,n)
    if (!is.na(first))
      v[1]=first
  }else{
    v[is.na(v)]=0.
  }
  v
}

#' Checks that the sum of specified fractions equals 1.0 
#' 
#' If all zero, replaces first element by 1.0
#' Generates an error if not all zero or sum != 1.0 
#' 
#' @param  v  a double[] vector
#' 
#' @return a vector with fractions that sum to 1.0 
check.sum.to.1 = function(v){
  if (abs(sum(v)-1.0)>1e-4)
    stop(paste("values of ",names(v)," should add up to 1.0",sep=""))
  v
}

#' Creates an approximate function to index a data frame using years specified in the rownames
#' 
#' @param  df      a data frame with row names that can be converted to numeric
#' @param  offset  an optional offset
#' 
#' 
#' @return an approxfun to generate a row index from a row name (usually a year)
create.approx.func.from.rownames=function(df,offset=0.0){
  x = as.numeric(rownames(df))+offset
  y = 1:nrow(df)
  if (nrow(df)==1){
    return(function(t){y})
  }else{
    return(approxfun(x,y,rule=2))
  }
}

#' Returns the deathrate for time t (read from WHO data using deathrate.from.data())
#' 
#' @param  p  params environment
#' @param  t  time
#' 
#' 
#' @return a vector of deathrates by age group at time t
get.deathrate.allcauses=function(p, t){
  return(p$deathmultiplier*p$deathrate[p$deathfn(t),])
}

#' Returns the background deathrate for time t based on a baseline run (see new.demography.from.model.run())
#' 
#' The background deathrate is the difference between the WHO deaths and the TB HIV deaths divided by population size.
#' See also the vignette on output.
#' 
#' @param  p  params environment
#' @param  t  time
#' 
#' 
#' @return a vector of background deathrate by age at time t
get.bgdeathrate=function(p, t){
  return(p$bgdeathrate[p$bgdeathfn(t),])
}

#' Creates an approximate function from an XML element and a data frame with parameter values
#' 
#' 
#' For example: 
#' <multiplier name="etamul" times="1960:2020" values="1/(1+exp(-sh*(x-midx)))"/>
#' 
#' @param  row a  named vector with names  'parameter', 'times', 'values' that will be parsed
#' @param  df  a  data frame with parameter values
#' 
#' 
#' @return an approxfun f(t)
create.approx.fun=function(row, df=NULL){
  s = as.character(row["values"]) ; values = unlist(strsplit(s,","))
  NAs = sum(is.na(suppressWarnings(as.numeric(values))))
  if (!is.null(df) & (NAs>0)){
    selcols  = which(names(df) %in% c("name","value"))
    df.left  = as.data.frame(df[,-selcols])
    sel = rowSums(is.na(df.left)) == ncol(df.left)
    df  = df[sel,selcols]
    eval(parse(text = paste0(df$name,"=",df$value)))
    s = unlist(strsplit(as.character(row["times"]),","))
    if (length(s)==1){
      eval(parse(text = paste0("x=",s)))
    }else{
      x = as.numeric(s)
    }
    y = 0 * x
    if (length(values)==1){
      eval(parse(text = paste0("y=",values)))
    }else{
      for (i in 1:length(y)){
        eval(parse(text = paste0("y[i]=",values[i])))
      }
    }
    y = as.numeric(y)
    return(approxfun(x,y,rule=2))
  }else{
  s = as.character(row["times"])
  x = as.numeric(unlist(strsplit(s,",")))
  y = as.numeric(values)
  return(approxfun(x,y,rule=2))
  }
}

#' Creates an approxfun of birthrate as a function of time from a filename and countrycode
#' 
#' @param  fname        filename
#' @param  countrycode  3 letter ISO country code (e.g. 'ZAF')
#' 
#' @return an approxfun birthrate(t)
calc.birthrate = function(fname=NA, countrycode=NA){
  if (!is.na(countrycode)){
    pop = read.csv(fname,header=F, fileEncoding = 'UTF-8-BOM')
    pop = pop[pop[,1]==countrycode,]
    nr  = nrow(pop)
    nc  = ncol(pop)
    x   = pop[2:nr,2]
    y   = pop[2:nr,3]/rowSums(pop[1:(nr-1),3:102])
    return(approxfun(x,y,rule=2))
  }
  return(NULL)
}

scale.population = function(pop=NA,fname=NA){
  if (!is.na(fname)){
    sumpop        = sum(pop)
    pop           = 0*pop
    pop.fracs     = as.matrix(read.delim(file=fname,sep="\t",header=T,stringsAsFactors = F, row.names=1, fileEncoding = 'UTF-8-BOM'))
    assert_that(dim(pop)[2]==dim(pop.fracs)[2])
    sum.pop.fracs = sum(pop.fracs)
    if (sum.pop.fracs<0.99 | sum.pop.fracs>1.01)
      warning(paste("The file ",fname, " does not appear to contain fractions adding up to 1.0 ?",sep=""))
    rowsel = rownames(pop) %in% rownames(pop.fracs)
    pop[rowsel,] = pop.fracs
    return(pop * sumpop)
  }
  NULL
}
#' Seeds the initial population with fractions by state as specified in the XML
#' 
#' @param  parms environment with parameters
#' 
#' @return an environment with parms$y.ini containing the initialized population
seed.initial.population=function(parms){
  z = t(parms$y.ini.read[1,] * t(parms$at.start))
  assert_that(sum(abs(colSums(z)-parms$y.ini.read[1,]))<1e-3,msg="seeded infections not specified correctly. NOTE: do not specify a fraction for the 1st stage !")
  parms$y.ini = z
}

#' Initializes constant parameters from XML
#' 
#' Called from read.model.parameters()
#' For example: 
#' <multiplier name="etamul" times="1960:2020" values="1/(1+exp(-sh*(x-midx)))"/>
#' 
#' @param  row a  named vector with names  'parameter', 'times', 'values' that will be parsed
#' @param  df  an approxfun f(t)  
#' 
#' 
#' @return an environment of model parameters to be passed to the next stage in parameter initialization
init.constants = function(p,xml){
  p$xml = xml
  with(p,{
    DAYSPERYEAR = 365.2425
    TB   = xml_attr(xml_find_all(xml$doc,"//TB/TB.stages/stage"),"name")     ; nTB   = length(TB) 
    HIV  = xml_attr(xml_find_all(xml$doc,"//HIV/HIV.stages/stage"),"name")   ; nHIV  = length(HIV) 
    RISK = xml_attr(xml_find_all(xml$doc,"//RISK/RISK.stages/stage"),"name") ; nRISK = length(RISK) 
    SES  = xml_attr(xml_find_all(xml$doc,"//SES/SES.stages/stage"),"name")   ; nSES  = length(SES) 
    VXa  = xml_attr(xml_find_all(xml$doc,"//VXa/VXa.stages/stage"),"name")   ; nVXa  = length(VXa)

    iTB   = check.sum.to.1(replace.NAs.by(xml_find_all(xml$doc,"//TB/TB.stages/stage"),"fraction.at.birth",0,1,nTB))  
    iHIV  = check.sum.to.1(replace.NAs.by(xml_find_all(xml$doc,"//HIV/HIV.stages/stage"),"fraction.at.birth",0,1,nHIV))     
    iRISK = check.sum.to.1(replace.NAs.by(xml_find_all(xml$doc,"//RISK/RISK.stages/stage"),"fraction.at.birth",0,1,nRISK))     
    iSES  = check.sum.to.1(replace.NAs.by(xml_find_all(xml$doc,"//SES/SES.stages/stage"),"fraction.at.birth",0,1,nSES))     
    iVXa  = check.sum.to.1(replace.NAs.by(xml_find_all(xml$doc,"//VXa/VXa.stages/stage"),"fraction.at.birth",0,1,nVXa))     
    
    ATBIRTH = list(VXa=iVXa,SES=iSES,RISK=iRISK,HIV=iHIV,TB=iTB)
    
    # NOTE : the order below VXa - SES - RISK - HIV - TB should NEVER be changed !!!
    
    DIMNAMES     = c("VXa","SES","RISK","HIV","TB")
    DIMNAMESLIST = list(VXa=VXa,SES=SES,RISK=RISK,HIV=HIV,TB=TB)
    DIMLENGTHS   = c(nVXa,nSES,nRISK,nHIV,nTB)
    names(DIMLENGTHS)=DIMNAMES
    
    HIVTB = NULL ; RISKHIVTB = NULL ; SESRISKHIVTB = NULL ; VXaSESRISKHIVTB = NULL;
    nHIV  = length(HIV)
    for (i in 1:nHIV) 
      HIVTB = c(HIVTB,paste(HIV[i],TB,sep="-")) 
    nHIVTB  = length(HIVTB)
    for (i in 1:nRISK) 
      RISKHIVTB = c(RISKHIVTB,paste(RISK[i],HIVTB,sep="-")) 
    nRISKHIVTB  = length(RISKHIVTB)
    for (i in 1:nSES)  
      SESRISKHIVTB = c(SESRISKHIVTB,paste(SES[i],RISKHIVTB,sep="-"))
    nSESRISKHIVTB  = length(SESRISKHIVTB)
    for (i in 1:nVXa)  
      VXaSESRISKHIVTB = c(VXaSESRISKHIVTB,paste(VXa[i],SESRISKHIVTB,sep="-"))
    nVXaSESRISKHIVTB  = length(VXaSESRISKHIVTB)

    inci.dim.names = list(VXa=calc.names.for.dim(environment(),"VXa"),
                         NULL,
                         RISK=calc.names.for.dim(environment(),"RISK"),
                         HIV=calc.names.for.dim(environment(),"HIV"),
                         TB=calc.names.for.dim(environment(),"TB"))

    eval(parse(text=paste("AGES=c(",xml_attr(xml_find_all(xml$doc,"//ages"),"lower.limits"),")",sep="")))
    names(AGES)      = paste("A",AGES,sep="") ; nAGES = length(AGES) 
    DEAD =  rep(F,nVXaSESRISKHIVTB)
    COUNT = rep(F,nVXaSESRISKHIVTB)
    for(name in names(p$DIMNAMESLIST)){
      DEAD  = DEAD  | grepl( "dead$",calc.names.for.dim(environment(),name))
      COUNT = COUNT | grepl("count$",calc.names.for.dim(environment(),name))
    }
    ALIVE     = !(DEAD | COUNT)   
  })
}

#' Creates contacts matrices 
#' 
#' See also the init.incidence() function in TBVx-data-reading-fncs-v##.R
#' 
#' @param  paths        paths to various files ; see set.paths() function 
#' @param  model.params params environment (to be updated)
#' 
#' @return a list(defined,M,MplustM) i.e. a logical to indicate if contacts matrices are defined and if so 
#' non NULL M and MplustM (i.e. M + t(M))
create.contacts.matrices=function(paths=NULL, model.params=NULL){
  assert_that(!is.null(paths),msg=paste("paths argument should not be NULL in",match.call()[[1]]))
  assert_that(is.environment(paths),msg=paste("paths argument should be an environment in",match.call()[[1]]))
  assert_that(!is.null(model.params),msg=paste("model.params argument should not be NULL in",match.call()[[1]]))
  if (!is.na(model.params$run.params$contact.matrix.txt)){
    temp=here(paths$country.dir,model.params$run.params$contact.matrix.txt)
    contactmatrixfname = stri_replace_last(str=temp,replacement=paste0("/",model.params$run.params$countrycode,"_"),fixed="/")
    M = expand.contact.matrix(model.params,read.contact.matrix(contactmatrixfname))
    MplustM = M+t(M)
    return(list(defined=T,M=M,MplustM=MplustM))
  }else{
    return(list(defined=F,M=NULL,MplustM=NULL))
  }
}

#' Initializes incidence data structures (TB, HIV and VXa) from XML and data files referenced in XML
#' 
#' See also the init.incidence() function in TBVx-data-reading-fncs-v##.R
#' 
#' @param  paths paths to various files ; see set.paths() function
#' @param  p     params environment to be updated
#' 
#' @return environment of model parameters with initialized incidence lists
initialize.incidence=function(paths=NULL,p=NULL){
  output = p$run.params$output
  inci.data = list()
  inci.data$TB   = inci.data.frame.from.xpath(p$xml,"//TB/TB.incidence/incidence.data")
  inci.data$HIV  = inci.data.frame.from.xpath(p$xml,"//HIV/HIV.incidence/incidence.data")      
  inci.data$RISK = inci.data.frame.from.xpath(p$xml,"//RISK/RISK.incidence/incidence.data")
  inci.data$VXa  = inci.data.frame.from.xpath(p$xml,"//VXa/VXa.incidence/incidence.data")
  
  for (dim in c("TB","HIV","RISK","VXa")){
    row = output$options$dim==dim
    tf  = output$options[row,"incidence"] & !is.null(inci.data[[dim]])
    output$options[row,"incidence"] = tf                                                
  }
  
  inci = list()
  for (name in names(inci.data)){
    df = inci.data[[name]]
    if (nrow(df)>0){
      inci[[name]]=list()
      for (i in seq_along(1:nrow(df))){
        inci[[name]][[i]]=init.incidence(p, fname=here(paths$country.dir,df$file[i]), 
                                            countrycode  = p$run.params$countrycode, 
                                            times=df$times[i], 
                                            values=df$values[i], 
                                            proportions=df$proportions[i], 
                                            denominator=df$denominator[i],
                                            once.per.year=df$once.per.year[i]==T,
                                            target.prevalence=df$target.prevalence[i]==T)
      }
    }
  }
  
  
  if (!is.na(p$intervention.start.from)){
    from = p$intervention.start.from
    rownr = 1
    assert_that(all(inci[[from]][[1]]$inci[rownr,]==0),msg=paste("the first row of the",from,"incidence data should only contain 0 values"))
    if (nrow(inci[[from]][[1]]$inci)>1){ rownr = 2 }
    p$intervention.start = as.numeric(rownames(inci[[from]][[1]]$inci)[rownr])
  }
  if (p$intervention){
    assert_that(p$intervention.start<2200,msg=paste("intervention start =",p$intervention.start,"; is that intended?"))
    modlog(level='WARN',paste0("intervention start = ",p$intervention.start," ??? ; is that intended?"))
  }
  inci
}

#' initializes TB transmission from the XML specification
#' 
#' @param  p   params environment
#' @return an environment to be assigned to Tm
initialize.TB.transmission=function(p=NULL){
  output = p$run.params$output
  xmltransitions = xml_find_all(p$xml$doc,"//TB/TB.transmission/transition.matrix/transition")
  if (length(xmltransitions)>0){    
    xmltimedparameters = xml_find_all(p$xml$doc,"//TB/TB.transmission/contact.rate.multiplier")
    xmlauxparameters   = xml_find_all(p$xml$doc,"//TB/TB.transmission/transition.matrix/parameter")
    
    Tm = init.parameters.from.xml(p,xmlpath="//TB/TB.transmission", xmlparameter="TB.parameter",
                                    xmltransitions=xmltransitions, xmlauxparameters=xmlauxparameters,xmltimedparams=xmltimedparameters)
    #Tm = optimize.params(Tm) 
    output$options[output$options$dim=="TB",]$transmission = T
    return(Tm)
  }else{
    return(NULL)
  }
}

#' initializes TB infectivity from the XML specification
#' 
#' @param  p   params environment
#' @return an environment to be assigned to Infc
initialize.TB.infectivity=function(p=NULL){
  xmltransitions = xml_find_all(p$xml$doc,"//TB/TB.infectivity/infectivity.matrix/infectivity")
  optimize.params(init.parameters.from.xml(p, xmlpath="//TB/TB.infectivity", xmlparameter="TB.parameter",xmltransitions=xmltransitions,values=T)) 
}

#' initializes HIV or TB treatment from the XML specification
#' 
#' @param  p   params environment
#' @param  dim a string to specify the dimension to initialize treatment for
#' @return an environment to be assigned to TBtr, HIVtr etc
initialize.treatment=function(p=NULL,dim=NA){
  if (p$DIMLENGTHS[dim]>1){
    output = p$run.params$output
    xmlpath = paste0("//",dim,"/",dim,".progression")
    nodeset = xml_find_all(p$xml$doc,paste0(xmlpath,"/treatment.matrix"))
    if (length(nodeset)>0){
      tr = list()
      for (i in 1:length(nodeset)){
        tr[[i]] = init.parameters.from.xml(p, xmlpath=xmlpath,
                                            xmlparameter=   paste0(dim,".parameter"),
                                            xmltransitions= xml_find_all(nodeset[i],"transition"),
                                            xmlauxparameters   = xml_find_all(nodeset[i],"parameter"),
                                            xmltimedparams= xml_find_all(nodeset[i],"multiplier"))
        names(tr)[i]=xml_attr(nodeset[i],"name")
        if (sum(is.na(as.matrix(tr[[i]]$parameter.matrices[[1]])))){
          print("found")
        }
        #tr[[i]] = optimize.params(tr[[i]]) 
        tr[[i]]$aggregated=F
      }
      output$options[output$options$dim==dim,]$treatment = T
      
      return(tr)
    }
  }
  return(NULL)
}

#' Initialize progression
#' 
#' initializes VXa, SES, RISK, HIV, or TB progression from the XML specification
#' 
#' @param  p   params environment
#' @param  dim (a string to specify the dimension to initialize progression for
#' @return an environment to be assigned to TBp, HIVp etc
initialize.progression=function(p=NULL,dim=NA){
  if (p$DIMLENGTHS[dim]>1){
    output = p$run.params$output
    xmlpath        = paste0("//",dim,"/",dim,".progression")
    xmltransitions = xml_find_all(p$xml$doc,paste0(xmlpath,"/transition.matrix/transition"))
    xmlauxparameters = xml_find_all(p$xml$doc,paste0(xmlpath,"/transition.matrix/parameter"))
    if (length(xmltransitions)>0){
      progression = init.parameters.from.xml(p,xmlpath=xmlpath,xmlparameter=paste0(dim,".parameter"),xmltransitions=xmltransitions, xmlauxparameters=xmlauxparameters)
      #progression = optimize.params(progression) 
      output$options[output$options$dim==dim,]$progression = T
      for (i in seq_along(p$inci[[dim]])){
        p$inci[[dim]][[i]]$inci.trend.fn = parse.time.series(p, p$inci[[dim]][[i]], progression$fixed.parameters)
      }
      progression$aggregated=F
      return(progression)
    }
  }
  return(NULL)
}

#' Reads model parameters from XML
#' 
#' reads the model parameters (without initialization ; called by the run() function) from XML
#' 
#' @param  paths an environment with paths to various files, country code etc. See the function set.paths()
#' @return an environment with initialized constants and data structures related to demography
#' @export
read.model.parameters=function(paths=NULL){
  assert <<- T
  assert_that(!is.null(paths),msg=paste("paths should not be NULL",match.call()[[1]]))
  assert_that(is.environment(paths),msg=paste("paths argument should be an environment in",match.call()[[1]]))
  assert_that(is.readable(paths$xml),msg=paste0("cannot read XML",paths$xml))

  p = new.env()
  with(p,{
    paths         = paths
    run.params    = parse.run.spec(paths$xml)
    if (!is.null(paths[["country.code"]])){ run.params$countrycode = paths$country.code }
    xmlfilename   = paths$xml
    init.constants(environment(),xml=read.xml(paths$xml,"/TB.Vx.model.inputfile"))
    intervention  = !is.na(paths[["baseline"]])
    aging.matrix  = create.aging.matrix(environment())
    if (!any(is.na(paths[["VXa.incidence.files"]]))){
      set.node.attrs(xml$doc, "//VXa/VXa.incidence/incidence.data", attrname="file", newvalues=paths$VXa.incidence.files)
    }
    assert_that(!is.na(run.params$population.csv),msg="path to population.csv is required")
    popfname = here(paths$country.dir,run.params$population.csv)
    y.ini.read = read.population(p           = environment(),
                            rownms      = p$VXaSESRISKHIVTB, 
                            colnms      = p$AGES, 
                            fname       = popfname, 
                            countrycode = run.params$countrycode, 
                            year        = run.params$simulation$years[1])
      
    assert_that(run.params$birthrate.from.data,msg="<birthrate from.population.data='true'/> is required in this model version")
    birthrate = calc.birthrate(fname=popfname, countrycode = run.params$countrycode)
    
    unpopdem       = get.demography(environment(),fname=popfname, countrycode=run.params$countrycode)
    #unpopdem       = apply(unpopdem,2,function(x){z = min(x[x>0]);x[x==0]=z;x})
    unpopdemfn     = create.approx.func.from.rownames(unpopdem)
    assert_that(!is.na(run.params$deathrate.txt),msg="path to deathrates.csv is required")
    deathrate      = deathrate.from.data(paths, environment())
    deathfn        = create.approx.func.from.rownames(deathrate)
    deathmultiplier = run.params$deathrate.multiplier
    contacts       = create.contacts.matrices(paths, environment())
    intervention.start.from = xml_attr(xml_find_all(xml$doc,"//simulation/options/intervention.start"),"from.incidence.data")
    if (length(intervention.start.from) == 0){
      intervention.start.from = NA
    }
    intervention.start = xml_attr(xml_find_all(xml$doc,"//simulation/options/intervention.start"),"year")
    if (length(intervention.start) == 0){
      intervention.start = 1e6
    }
    inci           = initialize.incidence(paths, environment()) # i.e. read the incidence data and initialize all but the (parameter dependent) values series
  })
  modlog(level = 'DEBUG',msg="successful result of read.model()")
  p
}

#' Sets population adjustment and background death rate in intervention equal to baseline run 
#' 
#' uses output$dPOPadj, output$popadjrateintv, output$dBGx, output$population to calculate
#' params$popadjrateintv, params$popadjfnintv, params$bgdeathrate and params$bgdeathfn
#' 
#' @param  params an environment with initialized parameters of the intervention run that will be modified
#' @param  output the output from the baseline run
#' 
#' @return NULL
set.baseline.in.model.parameters=function(params=NULL,output=NULL){
  assert_that(!is.null(output),msg="output argument should not be NULL in set.baseline.in.model.parameters()")
  assert_that(!params$intervention,msg="set baseline pop adjustment by EITHER providing a baseline arg pointing to files OR using the output data")
  params$popadjrateintv = matrix.from.popdf(params,pop=output$dPOPadj)
  params$popadjfnintv   = create.approx.func.from.rownames(params$popadjrateintv)
  Y                     = matrix.from.popdf(params,pop=output$dBGx)
  Z                     = matrix.from.popdf(params,pop=output$population)
  params$bgdeathrate    = Y/Z
  params$bgdeathfn      = create.approx.func.from.rownames(params$bgdeathrate)
  params$intervention=T
}

#' Combines TB progression, TB treatment, HIV progression, HIV treatment 
#' if matrices are constant over time and if age groups are the same.
#' Flags such as HIVp$aggregated and for each treatment matrix HIVtr[[i]]$aggregated 
#' keep track of aggregation status.
#' 
#' @param  TBp   TB progression matrices
#' @param  TBtr  TB treatment matrices
#' @param  HIVp  HIV progression matrices
#' @param  HIVtr HIV treatment matrices
#' 
#' @return aggregated TB and HIV progression and treatment matrices
initialize.aggregated.TBHIVptr=function(TBp=NULL,TBtr=NULL,HIVp=NULL,HIVtr=NULL){
  agg=new.env()
  agg$age.ranges = TBp$age.ranges
  agg$age.groups = TBp$age.groups
  agg$parameter.matrices = TBp$parameter.matrices
  TBp$aggregated=T
  for (i in seq_along(TBtr)){
    if (is.null(TBtr[[i]]$timed.parameters)){
      TBtr[[i]]$aggregated = T
      for (j in seq_along(TBtr[[i]]$parameter.matrices)){
        agg$parameter.matrices[[j]]=agg$parameter.matrices[[j]]+TBtr[[i]]$parameter.matrices[[j]]
      }
    }
  }
  if (is.null(HIVp)){
    return(agg)
  }
  if (all(TBp$age.groups==HIVp$age.groups)){
    HIVp$aggregated=T
    for (i in seq_along(HIVp$parameter.matrices)){
      agg$parameter.matrices[[i]]=agg$parameter.matrices[[i]]+HIVp$parameter.matrices[[i]]
    }
    for (i in seq_along(HIVtr)){
      if (is.null(HIVtr[[i]]$timed.parameters)){
        HIVtr[[i]]$aggregated = T
        for (j in seq_along(HIVtr[[i]]$parameter.matrices)){
          agg$parameter.matrices[[j]]=agg$parameter.matrices[[j]]+HIVtr[[i]]$parameter.matrices[[j]]
        }
      }
    }  
  }
  return(agg)
}

#' Initialization of the model parameters 
#' 
#' initializes the model parameters 
#' 
#' 
#' @param params An environment with the result of calling the function read.model.parameters(); 
#' the params argument contains e.g. data frames from files read but needs to be parsed further
#' @return an environment with initialized data structures required to run the simulation 
#' 
#' \strong{Description of the environment returned:}
#' \describe{
#' \item{\code{AGES: double [16] 0 5 10 15 ...}}{with lower limits of 5 year age groups}
#' \item{\code{aging.matrix: double [16 x 16] -0.2 0.2 0 0 0 0 0 0 0 0 ...   }}{a matrix used to move the population to the next higher age group once a year ; i.e. the result of multiplying the 16x16 aging matrix with the (in this case) 240 x 16 matrix of state variables}
#' \item{\code{ALIVE: logical [240] TRUE TRUE TRUE TRUE TRUE TRUE ...}}{logical vector that defines the ‘alive’ states, derived from state names not ending in ‘dead'}
#' \item{\code{at.birth: double [240] 0.4 0 0 0 0 0 0 0 0 0 ...}}{a vector with the distribution of newborns over states (sums to 1.0)}
#' \item{\code{ATBIRTH: list[5]}}{distributions of newborns by state used to calculate the \code{at.birth} vector ; for instance:}
#' \itemize{
#' \item \code{VXa : double [1] 1}
#' \item \code{SES : double [2] 0.4 0.6}
#' \item \code{RISK: double [1] 1}
#' \item \code{HIV : double [10] 1 0 0 0 0 0 0 0 0 0}
#' \item \code{TB  : double [12] 1 0 0 0 0 0 0 0 0 0 0 0}
#' }
#' \item{\code{at.start: double [240 x 16] 0.321475 0.037363 0.002855 0.030072 0.000481 ... ...}}{distribution of the population by state for each of the age groups at the start of the simulation ; each age group (column) sums to 1.0}
#' \item{\code{birthrate: function (v)}}{a function that calculates the row index from the current year}
#' \item{\code{contacts: list[3]}}{}
#' \itemize{
#' \item \code{defined: logical [1] TRUE}
#' \item \code{M      :Formal class 'dgeMatrix' [package "Matrix"] with 4 slots}
#' \item \code{MplustM:Formal class 'dgeMatrix' [package "Matrix"] with 4 slots}
#' \item M is the contacts matrix (from the paper by Kiesha Prem et al) and MplustM is M + t(M) which is used to rebalance the contacts matrix to the updated population composition which is done once per year. See the function update.contact.matrix(parms, alive.pop) in TBVx-age-related-fncs-v##.R)
#' }
#' \item{\code{DAYSPERYEAR: double [1] 365.2425}}{doubleber of days per year}
#' \item{\code{deathfn: function (v)}}{function to select the deathrate vector of the current year from the deathrate matrix (below)}
#' \item{\code{deathrate: double [151 x 16] 0.0527 0.0527 0.0527 0.0527 0.0512 ...}}{deathrate by year (row) and age group (column) derived from data (e.g. IND_deathrates.csv)}
#' \item{\code{DIMLENGTHS: integer [5] 1 2 1 10 12}}{lengths of the dimensions}
#' \item{\code{DIMNAMES: character [5] "VXa" "SES" "RISK" "HIV" "TB"}}{names of the dimensions}
#' \item{\code{DIMNAMESLIST: list[5]}}{list of state names by dimension}
#' \itemize{
#' \item \code{VXa : character [1] "never"}
#' \item \code{SES : character [2] "low" "high"}
#' \item \code{RISK: character [1] "risk1"}
#' \item \code{HIV : character [10] "HIV-" "HIVu1" "HIVu2" "HIVd1" ...}
#' \item \code{TB  : character [12] "Un" "Uc" "Lf" "Ls" ...}
#' }
#' \item{\code{hash: character [1] e.g. "9f1dfb8bec73a9aaa1a92064f9a924a9"}}{hash of this environment of initialized parameters}
#' \item{\code{HIV: character [10] "HIV-" "HIVu1" "HIVu2" "HIVd1" "HIVd2" "ARTn1" "ARTn2" "ARTs1" "ARTs2" ...}}{names of HIV states}
#' \item{\code{HIVp: environment [10]}}{environment with HIV progression parameters}
#' \item{\code{HIVTB: character [120] "HIV--Un" "HIV--Uc" "HIV--Lf" "HIV--Ls" "HIV--Ds" "HIV--Dc" "HIV--T" ...}}{names of the combination of HIV and TB states}
#' \item{\code{HIVtr: list[4]}}{list of HIV treatment environments; the names of these environments are derived from the XML}
#' \itemize{
#' \item \code{diag       : environment}
#' \item \code{treat      : environment}
#' \item \code{supp       : environment}
#' \item \code{HIVARTdead : environment}
#' }
#' \item{\code{inci: list[1]}}{list of incidence lists ; one for each supported dimension. Used for generating incidence matrices (function create.incidence.matrix()). Each of the list elements is an environment with all that is required to properly execute the transitions specified in the incidence files}
#' \item{\code{inci.dim.names: list[5]}}{list of dimensions that (currently) support incidence files}
#' \itemize{
#' \item \code{VXa: character [240] "never" "never" "never" "never" ...}
#' \item \code{[[2]]   : NULL}
#' \item \code{[[3]]   : NULL}
#' \item \code{HIV: character [240] "HIV-" "HIV-" "HIV-" "HIV-" ...}
#' \item \code{[[5]]   : NULL}
#' }
#' \item{\code{Infc: environment [11]}}{ Environment with parameters related to TB infectiousness ; most important is the matrix with infectiousness by state (a 240x240 matrix in this case)}
#' \item{\code{intervention: logical [1] FALSE}}{logical to flag intervention run}
#' \item{\code{intervention.start: double [1] 1e+06}}{default value of intervention start ; overwritten if applicable}
#' \item{\code{intervention.start.from: logical [1] NA}}{logical indicating if the intervention start is derived from incidence file}
#' \item{\code{iHIV: double [10] 1 0 0 0 0 0 0 0 0 0}}{distribution of newborns over HIV states}
#' \item{\code{iRISK :  double [1] 1}}{distribution of newborns over RISK states}
#' \item{\code{iSES :  double [2] 0.4 0.6}}{distribution of newborns over SES states}
#' \item{\code{iTB :  double [12] 1 0 0 0 0 0 0 0 0 0 ...}}{distribution of newborns over TB states}
#' \item{\code{iVXa :  double [1] 1}}{distribution of newborns over VXa states}
#' \item{\code{iHIV, iRISK, iSES, iTB, and iVXa }}{are used to generate ATBIRTH}
#' \item{\code{nAGES :  integer[1] 16}}{}
#' \item{\code{nHIV :  integer[1] 10}}{}
#' \item{\code{nHIVTB :  integer[1] 120}}{}
#' \item{\code{nRISK :  integer[1] 1}}{}
#' \item{\code{nRISKHIVTB :  integer[1] 120}}{}
#' \item{\code{nSES :  integer[1] 2}}{}
#' \item{\code{nSESRISKHIVTB :  integer[1] 240}}{}
#' \item{\code{nTB :  integer[1] 12}}{}
#' \item{\code{nVXa :  integer[1] 1}}{}
#' \item{\code{nVXaSESRISKHIVTB :  integer[1] 240}}{}
#' \item{}{nAGES, nHIV, nHIVTB, nRISK, nRISKHIVTB etceteta are the sizes of (combined) dimensions}
#' \item{\code{paths : environment [13]}}{environment with paths to various files ; see set.paths() function}
#' \item{\code{popfname :  character[1] e.g. "./countries/ZAF/data/demographics.csv"}}{filename of UNPOP demographics file (i.e. population development over time)}
#' \item{\code{PROGRESSION : list[5]}}{list with HIV, RISK, SES and VXa progression and in addition HIV treatment. Produced by PROGRESSION = list(HIV=HIVp,HIVtr=HIVtr,RISK=RISKp,SES=SESp,VXa=VXap)  where HIVp etc have been initialized separately. TB progression and treatment are in TBp and TBtr}
#' \itemize{
#' \item \code{HIV : environment}
#' \item \code{HIVtr: list[4]}
#' \item \code{RISK : NULL}
#' \item \code{SES  : NULL}
#' \item \code{VXa  : NULL}
#' }
#' \item{\code{RISK : character "risk1"}}{names of RISK states}
#' \item{\code{RISKHIVTB : character [120] "risk1-HIV--Un" "risk1-HIV--Uc" "risk1-HIV--Lf" "risk1-HIV--Ls" ...}}{names of combined RISKHIVTB states}
#' \item{\code{RISKp :  NULL}}{RISK progression}
#' \item{\code{run.params : environment}}{Environment with parameters that are related to running the simulation rather than the model itself (i.e. states, transitions etc.)}
#' \item{\code{SES :  character [2] "low" "high"}}{names of RISK states}
#' \item{\code{SESp :  NULL}}{SES progression}
#' \item{\code{SESRISKHIVTB :  character [240] "low-risk1-HIV--Un" "low-risk1-HIV--Uc" "low-risk1-HIV--Lf" ...}}{names of combined SESRISKHIVTB states}
#' \item{\code{TB :  character [12] "Un" "Uc" "Lf" "Ls" "Ds" "Dc" "T" "R" "TBdead" "TBHIVdead" "Rdead" ...}}{names of TB states}
#' \item{\code{TBHIVptr : environment[3]}}{TBHIV progression NOTE: to be explained!!}
#' \itemize{
#' \item \code{age.groups: integer[2]}
#' \item \code{age.ranges: list[2]}
#' \item \code{parameter.matrices : list[2]}
#' }
#' \item{\code{TBp : environment[10]}}{TB progression}
#' \item{\code{TBtr : list[5]}}{list of TB treatment matrices}
#' \itemize{
#' \item \code{RDs    : environment[10]} 
#' \item \code{LfLsDs : environment[10]}
#' \item \code{init   : environment[10]}
#' \item \code{nodead : environment[10]}
#' \item \code{died   : environment[10]}
#' }
#' \item{\code{Tm : environment[9]}}{TB transmission}
#' \item{\code{unpopdem :  double [151 x 16] 2024 2106 2218 2326 2419 ...}}{UNPOP population by year and age group (i.e. adjusted to modelled age groups)}
#' \item{\code{unpopdemfn : function (v)}}{function that returns the row index for the current year}
#' \item{\code{VXa :  character[1] "never"}}{names of VXa states}
#' \item{\code{VXap :  NULL}}{VXa progression}
#' \item{\code{VXaSESRISKHIVTB :  character [240] "never-low-risk1-HIV--Un" "never-low-risk1-HIV--Uc" ...}}{names of combined VXaSESRISKHIVTB states}
#' \item{\code{xml : list[2]}}{list with XML document and XML Schema}
#' \itemize{
#' \item \code{doc   :list[2]}
#' \item \code{schema:list[2]}
#' }
#' \item{\code{xmlfilename :  character[1] "./countries-examples/ZAF/parameters/XMLinput_2910_m.xml"}}{file name of XML input file}
#' \item{\code{y.ini :  double [240 x 16] 650.584 75.613 5.779 60.859 0.974 ...}}{initial population by state and age group}
#' \item{\code{y.ini.read :  double [240 x 16] 2024 0 0 0 0 ...}}{the initial population read from \code{demographics.csv} for the start year of the simulation, before seeding ; seeding creates \code{y.ini}}
#'}
#' @export
initialize.model.parameters=function(params=NULL){
   assert_that(!is.null(params),msg="pass (updated) parameters as an argument to initialize.model.parameters")
   with(params,{
    at.birth = fractions.by.stage.at.birth(environment())
    at.start = fractions.by.stage.at.start(environment(),df=get.seeded.infections(xml$doc))
    TBp   = initialize.progression(environment(),dim="TB")
    TBtr  = initialize.treatment(environment(),dim="TB")
    Tm    = initialize.TB.transmission(environment())
    Infc  = initialize.TB.infectivity(environment())
    HIVp  = initialize.progression(environment(),dim="HIV")
    HIVtr = initialize.treatment(environment(),dim="HIV")
    RISKp = initialize.progression(environment(),dim="RISK")
    RISKtr= initialize.treatment(environment(),dim="RISK")
    SESp  = initialize.progression(environment(),dim="SES")
    VXap  = initialize.progression(environment(),dim="VXa")
    TBHIVptr = initialize.aggregated.TBHIVptr(TBp=TBp,TBtr=TBtr,HIVp=HIVp,HIVtr=HIVtr)
    #TBHIVptr = optimize.params(TBHIVptr)
    seed.initial.population(environment())
  })
  params$hash  =  digest(params)
  modlog(level = 'DEBUG',msg="successful result of init.var.parameters.and.matrices()")
  params
}

