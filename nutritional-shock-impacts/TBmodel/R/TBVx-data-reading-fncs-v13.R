# Nov 14, 2022. Model version 3.2.8.9

#' Creates a file name based on the passed parameters, a hash of model parameters and the date and time
#' 
#' @param  path the path to which the file will be written
#' @param  model.params environment with initialized model parameters
#' @param  runnr   runnr (used when running the optimizer)
#' @param  runtype string which specifies the runtype e.g. "flows" or "stocks" etc
#' @param  score score (used when running the optimizer)
#' 
#' 
#' @return a file name that is usually unique due to including a hash of the model parameters and the date and time
#' @export
create.filename = function(path=NULL, model.params=NULL, runnr=NA, runtype=NA, score=NA){
  assert_that(!is.null(path),msg="path argument to create.filename should not be NULL")
  dattim  = format(Sys.time(), "%F_%Hh%Mm%Ss")
  s = paste0(path,"/[",model.params$run.params$countrycode,"]")
  if (exists(taskenvvar) & length(taskenvvar)>0){ 
    s = paste0(s,"[",taskenvvar,"]") 
  }
  s = paste0(s,"[",dattim,"]")
  s = paste0(s,"[",model.params$hash,"]") 
  charv = unlist(strsplit(model.params$xmlfilename,"/"))
  vlast = charv[length(charv)]
  s = paste0(s,"[",modelversion(),"][",unlist(strsplit(vlast,".xml"))[1],"]")
  if (!is.na(runnr)){ s = paste0(s,"[",runnr,"]")}
  if (!is.na(score)){ s = paste0(s,"[",score,"]")}
  if (!is.na(runtype)){s = paste0(s,"[",runtype,"]")}
  s
}


#' Reads and validates an XML file using the provided path to the XML Schema file
#' 
#' @param  xmlfile the name of the xml file
#' @param  schemapath the path to the XML schema (.xsd)
#' 
#' @return an XML document
#' @export
read.xml = function(xmlfile, schemapath){
  assert_that(is.readable(xmlfile),msg=paste0("reading the XML input file ", xmlfile," failed"))
  doc = read_xml(xmlfile)
  fpath=strsplit(xmlfile,"/",fixed=T)[[1]] # an array of the file path elements
  npath=paste(fpath[1:(length(fpath)-1)],collapse="/") # leave out the file name
  npath=paste0(npath,"/") # append a /
  xmlschemafile = paste0(npath,xml_attr(xml_find_first(doc,schemapath),"noNamespaceSchemaLocation"))
  assert_that(is.readable(xmlschemafile),msg=paste0("reading the XML Schema file ", xmlschemafile," failed"))
  assert_that(sum(is.na(stri_locate_last_fixed(xmlschemafile,schemaversion())))==0,msg=paste0("XML Schema file should be ",schemaversion()))
  schema = read_xml(xmlschemafile)
  valid = xml_validate(doc,schema)
  if (!valid){
    for (i in seq_along(attr(valid,"errors")))
      modlog(level='FATAL', msg= paste0("error when validating ",xmlfile," : ",attr(valid,"errors")[i]))
    stop(paste0("errors in ",xmlfile," ; see logfile for details"))
  }
  result = list(doc,schema)
  names(result)=c("doc","schema")
  result
}

#' Reads a tab delimited file with a contacts matrix 
#' 
#' @param  fname the name of a tab delimited file with 16 rows and 16 columns with names A0, A5, A10 ... A75
#' 
#' @return a 16 x 16 matrix 
#' @export
read.contact.matrix=function(fname){
  assert_that(is.readable(fname),msg="contacts matrix cannot be read")
  contacts=as.matrix(as.data.frame(read.table(fname,header=T)))
  assert_that(nrow(contacts)==ncol(contacts) & ncol(contacts)==16,msg="expecting a contacts matrix of 16 rows and 16 columns")
  assert_that(sum(colnames(contacts)==paste0("A",5*(0:15)))==16,msg=paste0("expecting column names to be ",paste0("A",5*(0:15))))
  rownames(contacts)=colnames(contacts)
  contacts
}

#' Reads a tab delimited file with population composition as produced by a previous simulation run (e.g. as in output$population)
#' 
#' 
#' Checks that age groups are the same as in the parameters environment
#' 
#' @param  fname the name of a tab delimited file with columns year, country code, and age group (0 thru 79, 80, 90) 
#' 
#' @return a matrix with population composition by year (in rownames) and age groups in columns
#' @export
matrix.from.popdf = function(p,pop=NULL){
  assert_that(!is.null(pop),msg="no pop df found in matrix.from.popdf")
  pop                    = as.data.frame(pop)
  sel                    = pop$country==p$run.params$countrycode 
  x                      = pop[sel,3:ncol(pop)]
  colnames(x)            = as.integer(colnames(x))
  nr                     = nrow(x)
  nk                     = ncol(x)
  assert_that(sum(colnames(x)==p$AGES)==nk,msg="columns in popdf not equal to age groups")
  rownames(x)            = pop$year[sel]
  return(as.matrix(x))
}

#' Reads a csv file with demographics data (from UNPOP) for a country by year (in rows) and single age years (in columns)
#' 
#' @param  p an environment with parameters that is used to convert the population data to the age groups in the parameters
#' @param  fname the file name
#' @param  countrycode the 3 letter ISO country code
#' 
#' @return a matrix with the population by year (rownames) and age group (colnames)
#' @export
get.demography = function(p,fname=NA, countrycode=NA){
  if (!is.na(countrycode) & !is.na(fname)){
    pop                    = read.csv(fname,header=F, fileEncoding = 'UTF-8-BOM')
    sel                    = pop[,1]==countrycode 
    x                      = pop[sel,3:102]
    colnames(x)            = c(0:99)
    nr                     = nrow(x)
    nk                     = ncol(x)
    y                      = aggregate.by.age.groups(x,p$AGES,sumcols=F,avg=F)
    colnames(y)            = names(p$AGES)
    rownames(y)            = pop[sel,][1:nr,2]    
    return(y)
  }
  return(NULL)
}

#' Reads a csv file with death rates(from UNPOP) for a country by year (in rows) and single age years (in columns)
#' 
#' @param  paths the path to which the country code is appended file name
#' @param  model.params an environment with parameters that is used to get the file name and country code from
#' 
#' @return a matrix with the death rates by year (rownames) and age group (colnames)
#' @export
deathrate.from.data = function(paths=NULL, model.params=NULL){
  
  assert_that(!is.null(paths),msg=paste("paths argument should not be NULL in",match.call()[[1]]))
  assert_that(is.environment(paths),msg=paste("paths argument should be an environment in",match.call()[[1]]))
  assert_that(!is.null(model.params),msg=paste("model.params argument should not be NULL in",match.call()[[1]]))
  
  temp        = here(paths$country.dir,model.params$run.params$deathrate.txt)
  countrycode = model.params$run.params$countrycode
  fname       = stri_replace_last(str=temp,replacement=paste0("/",countrycode,"_"),fixed="/")
  if (!is.na(countrycode) & !is.na(fname)){
    y                    = read.csv(fname,header=F, fileEncoding = 'UTF-8-BOM')
    sel                  = y[,1]==countrycode 
    x                    = y[sel,3:102]
    colnames(x)          = c(0:99)
    nr                   = nrow(x)
    nk                   = ncol(x)
    death.rate           = aggregate.by.age.groups(x,model.params$AGES,sumcols=F,avg=T)
    colnames(death.rate) = names(model.params$AGES)
    rownames(death.rate) = y[sel,][1:nr,2]    
    return(death.rate)
  }
  return(NULL)
}

#' Creates a matrix with an empty population 
#' 
#' @param  rownms rownames
#' @param  colnms colnames
#' 
#' @return a matrix with 0's in rows and columns 
#' @export
create.empty.population = function(rownms, colnms){
  pop = matrix(0.,nrow=length(rownms),ncol=length(colnms)) 
  colnames(pop)=colnms ; rownames(pop)=rownms
  pop
}

#' Reads the initial population composition from a file
#' 
#' @param  p      an environemnt with initialized parameters
#' @param  rownms rownames
#' @param  colnms colnames
#' @param  fname  the name of the .csv file
#' @param  countrycode 3 letter ISO country code
#' @param  year   the year (i.e. row in the file)
#' 
#' @return a matrix with the population composition in the year specified
#' @export
read.population = function(p, rownms, colnms, fname=NA, countrycode=NA, year=NA){
  if (!is.na(countrycode)){
    pop           = read.csv(fname,header=F, fileEncoding = 'UTF-8-BOM')
    sel           = pop[,1]==countrycode
    assert_that(nrow(pop[sel,])>=1,msg=paste0("no rows for country code ",countrycode," in ",fname))
    years         = pop[sel,2]
    assert_that(all(years == cummax(years)),msg=paste0("the 2nd column in ",fname," is supposed to contain monotonously increasing years ..."))
    if(ncol(pop) > 102){
      warning("ignoring columns from col nr 103 i.e. using columns 3 thru 102 for ages 0 thru 99")
    } 
    assert_that(ncol(pop)>=102,msg=paste("expecting at least 102 columns in ",fname," (ignoring columns from 103)"))
    z             = pop[sel,3:min(ncol(pop),102)]
    names(z)      = 0:99
    rownames(z)   = years
    z             = z[order(rownames(z)),]
    yrsel         = approxfun(rownames(z),1:nrow(z),rule = 2)    
    y             = aggregate.by.age.groups(z[yrsel(year),],p$AGES, sumcols = F, avg = F)
    x             = create.empty.population(rownms, colnms)
    x[1,]         = y  
    return(x)
  }
  return(NULL)
}

#' Reads a tab delimited incidence file 
#' 
#' @param  parms      an environemnt with initialized parameters
#' @param  fname  the name of the .txt file
#' @param  countrycode 3 letter ISO country code
#' @param  proportions a logical to indicate if the numbers in the file are proportions or rates (the default)
#' 
#' @return a list with dim, from, to, VXa, RISK, SES, HIV, TB, and inci ; see the vignette for more details on incidence files 
#' @export
read.incidence = function(parms, fname=NA, countrycode=NA, proportions=F){
  all.inci = read.delim(fname, header = T, stringsAsFactors = F, fileEncoding = 'UTF-8-BOM')
  country = unique(all.inci$country)
  assert_that(length(country)==1, msg = "values in 'country' column should be unique within file")
  assert_that(country==countrycode, msg = "'country' column should match country code")
  dim = unique(all.inci$dim)
  assert_that(length(dim)==1, msg = "values in 'dim' column should be unique within file")
  from = unique(all.inci$from)
  assert_that(length(from)==1, msg="values in 'from' column should be unique within file")
  to = unique(all.inci$to)
  assert_that(length(to)==1, msg="values in 'to' column should be unique within file")
  names(all.inci) = stri_replace_first_regex(names(all.inci),"^X","")
  sel.inci = all.inci
  rownames(sel.inci)=sel.inci$YEAR
  excl = names(sel.inci) %in% c("dim","from","to","country","YEAR","VXa","RISK","SES","HIV","TB")
  inci = as.matrix(sel.inci[,!excl])
  dimindex = which(parms$DIMNAMES==dim)
  
  compress = dim=="HIV" && 1==which(parms$DIMNAMESLIST[[dimindex]]==from)

  list(dim=dim,from=from,to=to,VXa =unique(all.inci$VXa),
                               RISK=unique(all.inci$RISK),
                               SES =unique(all.inci$SES),
                               HIV =unique(all.inci$HIV),
                               TB  =unique(all.inci$TB),
                               inci=adjust.age.groups(p=parms,M=inci,compress.1st=compress,proportions=proportions)) # compress incidence in newborns
}

#' Initializes an environment from an incidence file and options specified
#' 
#' 
#' The times and values arguments are converted to a function in the TBVx-initialization-funcs-v##.R
#' 
#' 
#' @param  p      an environment with initialized parameters
#' @param  fname  the name of the tab delimited .txt file with incidence 
#' @param  countrycode 3 letter ISO country code
#' @param  times  a series of times  
#' @param  values a series of values or an expression to be converted to a function producing values as a function of time 
#' @param  proportions a logical to indicate if the numbers in the file are proportions or rates (the default)
#' @param  denominator the denominator to be used in incidence calculations either 'susc' or 'all'
#' @param  once.per.year a logical to indicate if incidence numbers are rates (F) or fractions (T) that are moved to the new state once per year
#' 
#' @return an environment that is used in the derivs.inci.common() function in TBVx-derivs-v##.R with the inci.matrix as most important part
#' @export
init.incidence=function(p,fname=NA,countrycode=NA,times=NULL,values=NULL, proportions=F, denominator=NA, once.per.year=F, target.prevalence=F){
  e = new.env()
  with(e,{
    z           = read.incidence(p,fname,countrycode,proportions=proportions)
    dim         = z$dim
    assert_that(z$dim %in% names(p$inci.dim.names),msg=paste0("incidence data not supported for ",z$dim))
    froms       = stri_replace_all_regex(z$from,replacement = "",pattern="[:space:]")
    froms       = stri_split_fixed(froms,",")[[1]]
    check       = froms %in% p[[dim]]
    assert_that(all(check),msg = paste("not all from states: ",paste(froms,collapse=",")," found in ",p[[dim]]))
    tos         = stri_replace_all_regex(z$to,replacement = "",pattern="[:space:]")
    tos         = stri_split_fixed(tos,",")[[1]]
    check       = tos %in% p[[dim]]
    assert_that(all(check),msg = paste("not all to states: ",paste(tos,collapse=",")," found in ",p[[dim]]))
    # assert_that(from %in% p[[dim]],msg=paste(fname,":",from,"not in",dim,"stages defined in XML"))
    # assert_that(  to %in% p[[dim]],msg=paste(fname,":",  to,"not in",dim,"stages defined in XML"))
    assert_that((target.prevalence & once.per.year & !proportions & denominator=="susc") | (!target.prevalence),msg = paste("if target.prevalence == TRUE for ",p[[dim]],
                                                                                                                             "once.per.year and proportions should be F and denominator should be susc"))
    inci        = z$inci
    if (once.per.year){
      assert_that(sum(z$inci>1)==0,msg=paste(fname,": fraction in incidence file should be <=1 if once.per.year=true in XML"))
    }
    sel         = calc.names.for.dim(p,dim) %in% froms | calc.names.for.dim(p,dim) %in% tos
    susc        = calc.names.for.dim(p,dim) %in% froms & p$ALIVE # OK ; excludes COUNT as well
    tosel       = calc.names.for.dim(p,dim) %in% tos & p$ALIVE
    allsel      = p$ALIVE
    for (xdim in p$DIMNAMES){
      if (xdim != dim && !is.null(z[[xdim]]) && !is.na(z[[xdim]])){
        stages   = stri_trim_both(unlist(strsplit(z[[xdim]],",")))
        susc     = susc & calc.names.for.dim(p,xdim) %in% stages
        sel      = sel & calc.names.for.dim(p,xdim) %in% stages
        tosel    = tosel & calc.names.for.dim(p,xdim) %in% stages
        allsel   = allsel & calc.names.for.dim(p,xdim) %in% stages
      }
    }
    alive       = p$ALIVE
    inci.fn     = create.approx.func.from.rownames(inci,offset=1e-4)
    inci.matrix = create.incidence.matrix(p,dim,froms,tos,sel) 
    inci.times  = times
    inci.values = values
    inci.proportions = proportions
    inci.denominator = denominator
    inci.once.per.year = once.per.year
    inci.target.prevalence = target.prevalence
  })
  e
}