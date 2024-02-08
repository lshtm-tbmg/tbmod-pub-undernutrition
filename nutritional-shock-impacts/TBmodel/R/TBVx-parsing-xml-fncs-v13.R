# Nov 29, 2022. Model version 3.2.8.9									  

#' returns all possible attribute names for a specific XML parameter (e.g. "TB.parameter") from the XML Schema
#' 
#' @param xmlschema   the XML schema
#' @param xmlparameter the XML parameter
#' 
#' 
#' @return a character vector of attribute names e.g. "VXa.stage" "SES.stage" etc for xmlparameter "TB.parameter"
col.names.from.schema=function(xmlschema,xmlparameter){
  nodes=xml_find_all(xmlschema,"/xs:schema/xs:element")
  index=xml_attr(nodes,"name")==xmlparameter
  attributes = xml_find_all(nodes[index],"xs:complexType/xs:attribute")
  xml_attr(attributes,"name")
}

#' generates a query based on the column names passed and the index within each of the dimensions (VXa, SES, etc)
#' 
#' called by create.small.matrix.from.query()
#' 
#' @param parms     the initialized params environment 
#' @param col.names a vector of column names (e.g. "SES.stage","HIV.stage", etc) to generate the query from
#' @param vxa       index of the VXa stage to generate the query for
#' @param ses       index of the SES stage to generate the query for
#' @param risk      index of the RISK stage to generate the query for
#' @param hiv       index of the HIV stage to generate the query for
#' @param tb        index of the TB stage to generate the query for
#' @param ag        index of the age group to generate the query for
#' 
#' @return a query to be used in generating a small parameters matrix
query.from.col.names=function(parms=NULL, col.names=NULL, vxa, ses, risk, hiv, tb, ag){ # col.names should be derived from XML Schema !
  query=NULL
  for (name in col.names)
    query=paste(query,switch(
      name,
      "SES.stage" = paste("(is.na(SES.stage)  | SES.stage  == '",parms$SES[ses],"')",sep=""),
      "HIV.stage" = paste("(is.na(HIV.stage)  | HIV.stage  == '",parms$HIV[hiv],"')",sep=""),
      "TB.stage"  = paste("(is.na(TB.stage)   | TB.stage   == '",parms$TB[tb],"')",sep=""),
      "RISK.stage"= paste("(is.na(RISK.stage) | RISK.stage == '",parms$RISK[risk],"')",sep=""),
      "VXa.stage" = paste("(is.na(VXa.stage)  | VXa.stage  == '",parms$VXa[vxa],"')",sep=""),
      "age.group" = paste("(is.na(age.group)  | age.group  == '",ag,"')",sep="")),
      sep= " & ")
  substr(query,4,nchar(query))
}



#' parses an XML node list of <age.group name="A0"/> into a named integer vector of age.groups
#' 
#' @param xmldoc  the XML document
#' @param xmlpath the path to the XML element with <age.groups><age.group/>s
#' 
#' @return a named ('A0,'A15',..) integer vector (0,15,..)of age groups
parse.age.groupings.from.xml=function(xmldoc,xmlpath){
  nodes = xml_find_all(xmldoc,paste(xmlpath,"/age.groups/age.group"))
  if (length(nodes)==0){
    age.groups = 1
    names(age.groups)="A0"
    return(age.groups)
  }else{
    age.group.names=xml_attr(nodes,"name")
    assert_that(length(age.group.names)>0,msg=paste("no age group names found at ",xmlpath,sep=""))
    age.groups = as.integer(substring(age.group.names,2))
    names(age.groups)=age.group.names
    return(age.groups)
  }
}

#' verifies a data frame of fixed parameters for completeness of defined dependencies (i.e. all stages / age groups of a dimension)
#' 
#' @param p   the parameters environment
#' @param z   the data frame to be verified
#' @param xmlpath the XML path of the parameters
#' @param age.groups the age groups defined 
#' 
#' @return NULL
verify.data.frame.with.fixed.parameters=function(p, z, xmlpath, age.groups=NULL){
  parameters = unique(z$name)
  y = z[,-which(names(z)=="value")]
  for (parameter in parameters){
    x = y[y$name==parameter,]
    nrows = nrow(x)
    for (colname in names(x)[-which(names(x)=="name")]){
      nNAs = sum(is.na(x[,colname]))
      assert_that(!(nNAs>0 & nNAs<nrows),msg=paste("error at ",xmlpath," in specification of parameter ",parameter,": specify dependency on ",colname," either as NA or specific for all ",colname,sep=""))
      if (nNAs==0){
        if (colname!="age.group"){
          dim = stri_split_fixed(colname,".")[[1]][1]
          assert_that(all(unique(x[,colname]) %in% p[[dim]]),msg=paste("specify value of",parameter,"for all",dim,"stages"))
          assert_that(0==(length(x[,colname]) %% p$DIMLENGTHS[[dim]]),msg=paste("specify value of",parameter,"once for all",dim,"stages"))
        }else{
          assert_that(all(unique(x[,colname]) %in% names(age.groups)),msg=paste("specify value of",parameter,"for all age groups"))
          assert_that(0==(length(x[,colname]) %% length(age.groups)),msg=paste("specify value of",parameter,"once for all age.groups"))
        }
      }
    }
  }
}

#' parses XML timed parameters such as <multiplier name="etamul" times="1960:2020" values="1/(1+exp(-sh*(x-midx)))"/> <multiplier>
#' 
#' @param p   the parameters environment
#' @param xmlnodes a node list of XML parameters with name, times and values attributes
#' @param fixed.params fixed parameters used in parsing the 'values' expression
#' 
#' @return an environment with named functions of time 
parse.timed.parameters=function(p,xmlnodes, fixed.params){
  df = data.frame(parameter=xml_attr(xmlnodes,"name"),times=xml_attr(xmlnodes,"times"),values=xml_attr(xmlnodes,"values"),stringsAsFactors = F)
  if (nrow(df)>0){
    if (nrow(df)>1 & sum(duplicated(df$parameter))>0)
      stop(paste("duplicate definition of time dependent parameters ; please correct input file ... exiting"))
    fns=apply(df,1,create.approx.fun, fixed.params)
    e = new.env()
    for (i in seq_along(df$parameter))
      eval(parse(text=paste(df$parameter[i],"=fns[[i]]",sep="")),envir=e)
    return(e)
  }
}

#' parses a series of times and values derived from the incidence list elements inci.values and inci.times 
#' 
#' @param p   the parameters environment
#' @param inci a list of incidence parameters
#' @param fixed.params fixed parameters used in parsing the 'values' expression
#' 
#' @return a function of time for multiplying the incidence with
parse.time.series=function(p,inci, fixed.params){
  arg = list(times=inci$inci.times,values=inci$inci.values)
  fnc = create.approx.fun(arg,fixed.params)
  return(fnc)
}

#' checks that all XML parameters at a specified xmlpath have a name and a value attribute
#' 
#' @param xmlpath   the path to the XML parameters
#' @param xmlparams a node list of XML parameters
#' 
#' @return NULL
check.names.values=function(xmlpath,xmlparams){
  n = length(xmlparams)
  names = xml_attr(xmlparams,"name")
  values = xml_attr(xmlparams,"value")
  assert_that(n>0,msg=paste("no parameters found at ",xmlpath,sep=""))
  assert_that(length(names)==length(values), msg=paste("inconsistent or incomplete parameters at ",xmlpath,sep=""))
}

#' parse fixed parameter values from XML into a data frame
#' 
#' @param p initialized parameters environment
#' @param xmlpath path to the XML elements e.g. "//TB/progression"
#' @param xmlparameter name of the parameter e.g. "TB.parameter"
#' @param col.names  column names to use for parameter value dependencies
#' @param age.groups a named vector of age group lower bounds
#' @param dim.name the dimension e.g. "TB"
#' 
#' @return a data frame with parameter names and values including dependencies ; column names as passed in col.names
parse.fixed.parameters.from.xml=function(p,xmlpath,xmlparameter,col.names=NULL,age.groups=NULL,dim.name){
  if (length(xml_find_all(p$xml$doc,xmlpath))<=0)
    return(NULL)
  xmlparams  = xml_find_all(p$xml$doc,paste(xmlpath,xmlparameter,sep="/"))
  check.names.values(xmlpath,xmlparams)
  X = list()
  for (i in seq_along(col.names))
    X[[i]]=xml_attr(xmlparams,col.names[i])
  names(X)=col.names
  df = data.frame(X,stringsAsFactors = F)
  verify.data.frame.with.fixed.parameters(p, df, xmlpath, age.groups)
  df  
}

#' creates a small parameter matrix for an age group, a dimension (e.g. "TB") and a selection of index values of the other dimensions from a query

#' @param parms_ initialized parameters environment
#' @param query_ the query generated by query.from.col.names()
#' @param z_ a data frame of parameter values (by name and value) including dependencies (for all dimensions and age groups except current dimension)
#' from which a subset will be generated using the query
#' @param m_ a data frame with transition rates (if values=F) or values by stage (if values=T e.g. for infectivity)
#' @param auxparameters_ auxiliary parameters defined in the transition or treatment matrix (usually not used)
#' @param values_ a logical to indicate if values should be read for stages (e.g. for TB infectivity) or for transitions (from, to)
#' @param dim.name_ the dimension e.g. "TB" 
#' 
#' @return a double p$DIMLENGTHS x p$DIMLENGTHS matrix with parameter values to be used to calculate the contribution of e.g. progression to dY/dt
create.small.matrix.from.query=function(parms_, query_=NA, z_, m_, auxparameters_, values_, dim.name_){
  index_ = which(parms_$DIMNAMES==dim.name_)
  M_     = matrix(0,parms_$DIMLENGTHS[index_],parms_$DIMLENGTHS[index_])
  colnames(M_)=rownames(M_)=parms_$DIMNAMESLIST[[index_]]
  if (!is.na(query_) & nrow(z_)>0){
    z_=z_[with(z_,eval(parse(text=query_))),]
  }
  if (nrow(z_)>0) # no need to iterate over rows of data frame z :-)
    eval(parse(text=paste(z_$name,"=",z_$value)))
  if (!is.null(auxparameters_) && nrow(auxparameters_)>0) # no need to iterate over rows of data frame z :-)
    eval(parse(text=paste(auxparameters_$name,"=",auxparameters_$value)))
  for (row_ in 1:nrow(m_)){
    if (!values_){
      #print(paste("M_=transition(M_,",'m_$from[row_]',",",'m_$to[row_]',",",m_$rate[row_],")",sep=""))
      eval(parse(text = paste("M_=transition(M_,",'m_$from[row_]',",",'m_$to[row_]',",",m_$rate[row_],")",sep="")))
    }else{
      eval(parse(text=paste("M_[",'m_$stage[row_]',",",'m_$stage[row_]',"]=",m_$value[row_],sep="")))
    }
  }
  M_
}

#' creates a parameter matrix for an age group e.g. for progression, infectivity, etc
#' 
#' 
#' calls  create.small.matrix.from.query()
#' 

#' @param p initialized parameters environment
#' @param ag age group 
#' @param params.df a data frame of parameter values (by name and value) including dependencies (for all dimensions and age groups except current dimension)
#' @param transitions a data frame of transition with columns from, to and rate
#' @param auxparameters auxiliary parameters defined in the transition or treatment matrix (usually not used)
#' @param values a logical to indicate if values should be read for stages (e.g. for TB infectivity) or for transitions (from, to)
#' @param dim.name the dimension e.g. "TB" 
#' 
#' @return a double p$DIMLENGTHS x p$DIMLENGTHS matrix with parameter values to be used to calculate the contribution of e.g. progression to dY/dt
create.full.parameter.matrix.for.age.group = function(p=NULL,ag=NULL,params.df,transitions,auxparameters,values,dim.name){
  sub.df = NULL
  if (nrow(params.df)<2){
    sub.df = params.df[,apply(params.df,2,is.na)!=nrow(params.df)]
  }else{
    sub.df = params.df[,apply(apply(params.df,2,is.na),2,sum)!=nrow(params.df)]
  }  
  largeM = matrix(0,prod(p$DIMLENGTHS),prod(p$DIMLENGTHS))
  smallM = NULL ; query  = NULL
  nTB    = p$nTB  ; nHIV = p$nHIV ; nSES = p$nSES ; nRISK = p$nRISK ; nVXa = p$nVXa
  switch (dim.name, "HIV" = {nHIV = 1}, "TB" = {nTB = 1}, "VXa" = {nVXa = 1}, "SES" = {nSES = 1}, "RISK" = {nRISK = 1})
  col.names = names(sub.df)[-c(which(names(sub.df)=="name"),which(names(sub.df)=="value"))]
  
  for (vxa in 1:nVXa){
    for (hiv in 1:nHIV){ 
      for (risk in 1:nRISK){ 
        for (ses in 1:nSES){ 
          for (tb in 1:nTB){
            sel = as.list(c(vxa,ses,risk,hiv,tb))
            names(sel)=p$DIMNAMES
            sel[[dim.name]]=1:p$DIMLENGTHS[dim.name]
            i = calc.small.matrix.indices.for.one.full.dim.and.specific.single.dim.indices.from.list(p,sel)
            if (length(col.names)==0){
              if (is.null(smallM))
                smallM = create.small.matrix.from.query(p, NA, sub.df, transitions, auxparameters, values, dim.name)
              largeM[i,i] = smallM
            }else{
              newquery = query.from.col.names(p, col.names, vxa, ses, risk, hiv, tb, ag)
              if (is.null(query) || stri_cmp(query,newquery)!=0){
                smallM = create.small.matrix.from.query(p, newquery, sub.df, transitions, auxparameters, values, dim.name)
              }
              largeM[i,i] = smallM
              query = newquery
            }
          }
        }
      }
    }
  }
  largeM[,!p$ALIVE]=0. # you can never leave the dead state .... at least not by transition
  Matrix(largeM)
}

#' creates a list of parameter matrices by age group e.g. for progression, infectivity, etc
#' 
#' 
#' calls  create.full.parameter.matrix.for.age.group()
#' 
#' @param p initialized parameters environment
#' @param age.groups 
#' @param fixed.parameters an environment with xmlpath, xmlparameter, dim.name, col.names, age.groups, age.ranges, parameter.matrices, timed.parameters
#' @param xmlpath a node list of XML parameters defined in the transition matrix (if any) i.e. "//TB/TB.transmission/transition.matrix/parameter"
#' @param xmltransitions a node list of XML timed parameters i.e. time dependent multipliers from e.g. "//TB/TB.transmission/transition.matrix/multiplier"
#' @param xmlauxparameters a logical to indicate if values should be read for stages (e.g. for TB infectivity) instead of for transitions (from, to)
#' @param values a logical to indicate if values should be read for stages (e.g. for TB infectivity) or for transitions (from, to)
#' @param col.names the column names of the transition matrix
#' @param dim.name the dimension e.g. "TB" 
#' 
#' @return a list (by age group) of parameter matrices for the specified dimension
create.list.of.parameter.matrices.by.age.group = function(p,age.groups,fixed.parameters,xmlpath,xmltransitions,xmlauxparameters,values,col.names=NULL,dim.name=NULL){
  n  = length(age.groups)
  assert_that(n>0)
  L  = list()
  transitions = NULL
  auxparameters = NULL
  if (!values){
    nodes=xmltransitions
    assert_that(length(nodes)>0)
    transitions = data.frame(from=xml_attr(nodes,"from"),to=xml_attr(nodes,"to"),rate=xml_attr(nodes,"rate"),stringsAsFactors = F)
    nodes=xmlauxparameters
    if (!is.null(nodes) && length(nodes)>0)
      auxparameters = data.frame(name=xml_attr(nodes,"name"),value=xml_attr(nodes,"value"),stringsAsFactors = F)
  }else{
    nodes=xmltransitions
    assert_that(length(nodes)>0)
    transitions = data.frame(stage=xml_attr(nodes,"stage"),value=xml_attr(nodes,"value"),stringsAsFactors = F)
  }
  for (a in 1:n){
    L[[a]] = create.full.parameter.matrix.for.age.group(p,ag=names(age.groups)[a],fixed.parameters,transitions,auxparameters,values,dim.name)
    if (sum(is.na(L[[a]]))>0){
      print("gotcha")
    }
  }
  L  
}
#' initializes an environment with parameters for infectivity, transmission (Tm), progression, treatment
#' 
#' @param xmlpath   the path to the XML parameters e.g. "//TB/TB.transmission"
#' @param xmlparameter the name of the XML parameter elements e.g. "TB.parameter" 
#' @param xmltransitions a node list of XML transitions e.g. from "//TB/TB.transmission/transition.matrix/transition"
#' @param xmlauxparameters a node list of XML parameters defined in the transition matrix (if any) i.e. "//TB/TB.transmission/transition.matrix/parameter"
#' @param xmltimedparameters a node list of XML timed parameters i.e. time dependent multipliers from e.g. "//TB/TB.transmission/transition.matrix/multiplier"
#' @param values a logical to indicate if values should be read for stages (e.g. for TB infectivity) instead of for transitions (from, to)
#' 
#' @return an environment with xmlpath, xmlparameter, dim.name, col.names, age.groups, age.ranges, parameter.matrices, timed.paramters
init.parameters.from.xml=function(p,xmlpath,xmlparameter,xmltransitions,xmlauxparameters,xmltimedparams=NA,values=F){
  e = new.env()
  with(e,{
    xmlpath            = xmlpath
    xmlparameter       = xmlparameter
    dim.name           = strsplit(xmlparameter,".",fixed=T)[[1]][1]
    col.names          = col.names.from.schema(p$xml$schema,xmlparameter)
    age.groups         = parse.age.groupings.from.xml(p$xml$doc,xmlpath)
    age.ranges         = get.age.indices(p,age.groups)
    fixed.parameters   = parse.fixed.parameters.from.xml(p,xmlpath,xmlparameter,col.names,age.groups,dim.name)
    parameter.matrices = create.list.of.parameter.matrices.by.age.group(p,age.groups,fixed.parameters,xmlpath,xmltransitions,xmlauxparameters,values,col.names,dim.name)
    if (!is.na(xmltimedparams) && length(xmltimedparams)>0){
      timed.parameters = parse.timed.parameters(p,xmltimedparams, fixed.parameters)
    }else{
      timed.parameters = NULL
    }
  })
  e  
}

# not used
#default.value.from.schema.by.element.and.attribute=function(xmlschema,elementname,attrname,type){
#  j = which(xml_attr(xml_find_all(xmlschema,"//xs:element"),"name")==elementname)
#  element = xml_find_all(xmlschema,"//xs:element")[j]
#  k = which(xml_attr(xml_find_all(element,"xs:complexType//xs:attribute"),"name") == attrname)
#  default = xml_attr(xml_find_all(element,"xs:complexType//xs:attribute")[k],"default")
#  if (type=="double"){
#    return(as.double(default))
#  }else if (type=="integer"){
#    return(as.integer(default))
#  }else{
#    return(default)
#  }
#  return(NULL)
#}

#' retrieves the default value of an XML attribute from the XML Schema
#' 
#' @param xmlschema   the XML schema
#' @param attrname the attribute name
#' @param type the type of the attribute (e.g. double, fraction, integer, ...)
#' 
#' @return the default value of the attribute passed as an argument
default.value.from.schema=function(xmlschema,attrname,type){
  named.elements = xml_attr(xml_find_all(xmlschema,"//xs:element"),"name")
  sel = xml_attr(xml_find_all(xmlschema,"//xs:element//xs:attribute"),"name") == attrname
  default = xml_attr(xml_find_all(xmlschema,"//xs:element//xs:attribute")[sel],"default")
  if (type=="double" | type=="fraction"){
    return(as.double(default))
  }else if (type=="integer"){
    return(as.integer(default))
  }else{
    return(default)
  }
  return(NULL)
}

#' reads the attributes from the XML element <simulation><options><numerical.integration/> including defaults from the XML Schema
#' to be stored in params$num.int 
#' 
#' @param xmlschema   the XML schema
#' @param node a nodelist
#' @param names a vector of attribute names
#' @param types a vector of parameter types
#'  
#' 
#' @return a named list of values about the numerical integration algorithm
rk.parameters = function(xmlschema,node=NA,names=NA,types=NA){
  p = list()
  for (i in seq_along(names)){
    nm  = names[i]
    typ = types[i]
    q = xml_attr(node,nm)
    if (is.na(q)){
      p[[nm]] = default.value.from.schema(xmlschema,nm,typ)
    }else{
      if (typ=="double" | typ=="fraction"){
        p[[nm]] = as.double(q)
      }else if (typ == "integer"){
        p[[nm]] = as.integer(q)
      }else{
        p[[nm]] = q
      }
    }
  }
  p
}

#' reads the attributes from the XML element <simulation>
#' to be stored in params$simulation
#' 
#' @param xmldoc the XML document
#' @param h dt
#' 
#' @return a named list of start, stop, and years (intermediate output moments) all in years
simulation.years= function(xmldoc=NULL, h=NULL){
  simulation = xml_find_all(xmldoc,"//simulation")
  start = as.integer(xml_attr(simulation,"from.year"))
  stop  = as.integer(xml_attr(simulation,"thru.year"))
  list(start=start,stop=stop,years=start+cumsum(c(0,rep(h,(stop-start)/h))))
}

#' reads the XML file and parses parameters unrelated to the model itself (i.e. dimensions, progression etc etc)
#' 
#' @param xmlfile the XML file
#' 
#' @return an environment the will be part of the initialized parameters environment as run.params
parse.run.spec = function(xmlfile=NA){
  
  p = new.env()
  p$population.csv = NA
  p$countrycode = NA
  p$birthrate.from.data = F
  p$birthrate = NA
  p$birthrate.is.number = F
  p$deathrate.from.file = F
  p$deathrate.txt = NA
  p$contact.matrix.txt = NA
  p$population.fractions = NA
  p$population.total = NA
  
  xmllist              = read.xml(xmlfile,"/TB.Vx.model.inputfile")
  xmldoc               = xmllist$doc
  xmlschema            = xmllist$schema

  p$countrycode        = xml_attr(xml_find_all(xmldoc,"//simulation/demography"),"country.code")
  p$rescale.pop.yr     = as.numeric(xml_attr(xml_find_all(xmldoc,"//simulation/demography"),"rescale.population"))
  p$rescaled.pop       = F
  demography.from.data = xml_find_all(xmldoc,"//simulation/demography/from.data")
  if (length(demography.from.data)>0){
    xmlpopulation     = xml_find_all(demography.from.data,"//population")
    p$population.csv  = xml_attr(xmlpopulation,"file")
    xmlbirthrate      = xml_find_all(demography.from.data,"//birthrate")
    p$birthrate.from.data = xml_attr(xmlbirthrate,"from.population.data") == "true"
    if (!p$birthrate.from.data){
      p = parse.birthrate(p,xmlbirthrate)
    }
    xmlmortality      = xml_find_all(demography.from.data,"//mortality")
    if (!is.na(xml_attr(xmlmortality,"file"))){
      p$deathrate.from.file = T
      p$deathrate.txt = xml_attr(xmlmortality,"file")
      p$deathrate.multiplier = as.numeric(ifelse(is.na(xml_attr(xmlmortality,"multiplier")),1.0,xml_attr(xmlmortality,"multiplier")))
    }
  }
  
  demography.from.dist = xml_find_all(xmldoc,"//simulation/demography/from.distribution")
  if (length(demography.from.dist)>0){
    xmlpopulation          = xml_find_all(demography.from.dist,"//population")
    p$population.fractions = xml_attr(xmlpopulation,"file")
    p$population.total     = as.numeric(xml_attr(xmlpopulation,"total"))
    p                      = parse.birthrate(p,xml_find_all(demography.from.dist,"//birthrate"))
    p$deathrate.txt        = xml_attr(xml_find_all(demography.from.dist,"//mortality"),"file")
  }
  
  xml.contact.matrix   = xml_find_all(xmldoc,"//simulation/contact.matrix")
  p$contact.matrix.txt = ifelse(length(xml.contact.matrix)>0,xml_attr(xml.contact.matrix,"file"),NA)

  p$incidence.files    = as.list(xml_attr(xml_find_all(xmldoc,"//simulation/incidence/data"),"file"))
  

  p$dt = as.double(xml_attr(xml_find_all(xmldoc,"//simulation"),"dt"))
  if (is.na(p$dt)){
    p$dt = default.value.from.schema(xmlschema,"dt","fraction")
  }
  nodes = xml_find_all(xmldoc,"//simulation/options/numerical.integration")
  for (i in seq_along(nodes)){
    p$num.int[[i]]=rk.parameters(xmlschema,node=nodes[i],
                                 names=c("from.year","method","maxsteps","hmin","hini","atol","rtol","min.value.for.state.var"),
                                 types=c("double","text","integer","double","double","double","double","double"))
  }  

  finpopfrac = xml_find_all(xmldoc,"//output/final.population.as.fraction")
  if (length(finpopfrac)>0 && xml_attr(finpopfrac,"file")!="NA"){
    p$final.population.as.fraction = xml_attr(finpopfrac,"file")
  }else{
    p$final.population.as.fraction = NA
  }
  p$simulation = simulation.years(xmldoc,p$dt)
  p$output     = parse.output.spec(xmldoc)
  assert_that(all(p$output$years %in% p$simulation$years),msg=" not all years specified for output within simulation years")

  p
}

#' reads the fractions of seeded TB infections by stage from <TB><seeded.infections>
#' 
#' @param xmldoc the XML document
#' 
#' @return a data frame with seeded infections ; columns stage, age.from, age.thru, fraction
get.seeded.infections = function(xmldoc){
  seeded.infections = xml_find_all(xmldoc,"//TB/seeded.infections/seeded.infection")
  col.names = c("stage","age.from","age.thru","fraction")
  X = list()
  X[["stage"]]=xml_attr(seeded.infections,"stage")
  X[["age.from"]]=as.integer(xml_attr(seeded.infections,"age.from"))
  X[["age.thru"]]=as.integer(xml_attr(seeded.infections,"age.thru"))
  X[["fraction"]]=as.numeric(xml_attr(seeded.infections,"fraction"))
  df = data.frame(X,stringsAsFactors = F)
  multiplier = as.numeric(xml_attr(xml_find_all(xmldoc,"//TB/seeded.infections"),"value"))
  fromages = unique(df$age.from)
  thruages = unique(df$age.thru)
  uninfected.stage = xml_attr(xml_find_all(xmldoc,"//TB/TB.stages/stage")[1],"name")
  z  = NULL
  for (i in seq_along(fromages)){
    sel  = df$age.from == fromages[i]
    df$fraction[sel] = df$fraction[sel] * multiplier
    assert_that(sum(df$fraction[sel])<=1,msg="multiplier too large, sum of initial TB stage fractions > 1")
    assert_that(unique(df$age.thru[sel])==thruages[i],msg="inconsistency in from and thru ages in seed.infections")
    fraction.uninfected  = 1 - sum(df$fraction[sel])  
    z  = rbind(z,data.frame(stage=uninfected.stage,age.from=fromages[i],age.thru=thruages[i],fraction=fraction.uninfected))
  }
  rbind(df,z)  
}

# not used
# parse.birthrate=function(p,xmlbirthrate){
#  cbr.from.pop.data = xml_attr(xmlbirthrate,"from.population.data")
#  if (!is.na(cbr.from.pop.data) & cbr.from.pop.data!="false"){
#    stop("birthrate attribute 'from.population.data' should be 'false' when parsing other options")
#  }
#  cbr.file     = xml_attr(xmlbirthrate,"file")
#  cbr.number   = xml_attr(xmlbirthrate,"number")
#  cbr.fraction = xml_attr(xmlbirthrate,"fraction")
#  if (sum(c(is.na(cbr.fraction),is.na(cbr.file),is.na(cbr.number)))!=2)
#    stop("specify either a file, a fraction or a number for crude birth rate")
#  
#  p$birthrate.is.number=F
#  
#  if (!is.na(cbr.file)){
#    df = read.delim(cbr.file,stringsAsFactors = F, fileEncoding = 'UTF-8-BOM')
#    if (names(df)!=c("YEAR","CBR")){
#      if (min(df[,1])>=1800 & max(df[,1])<=2200 & min(df[,2])>=0. & max(df[,2])<0.2){
#        print("assuming 1st column of birthrate data contains the year and 2nd column contains the crude birth rate as a fraction")
#        p$birthrate = approxfun(x=df[,1],y=df[,2],rule=2)
#      }else{
#        stop("cannot parse birth rate data ; please provide data in a tab-delimited file with columns named YEAR and CBR")
#      }
#    }else{
#      p$birthrate = approxfun(x=df$YEAR,y=df$CBR,rule=2)
#    }
#  }
#  if (!is.na(cbr.fraction)){
#    p$birthrate = approxfun(x=c(1800,2200),y=rep(as.numeric(cbr.fraction),2),rule=2)
#  }
#  if (!is.na(cbr.number)){
#    p$birthrate.is.number = T
#    p$birthrate = approxfun(x=c(1800,2200),y=rep(as.numeric(cbr.number),2),rule=2)
#  }
#  p     
#}

#' parses the contents of the XML element <output><detailed.output> 
#' 
#' @param xmldoc the XML document
#' @param xmlpath by default "//output/detailed.output"
#' 
#' @return an environment with details on which output to generate and when ; see XML
parse.output.spec=function(xmldoc,xmlpath="//output/detailed.output"){
  if (length(xml_find_all(xmldoc,xmlpath))<=0){
    return(NULL)
  }  
  e = new.env()
  detailed   = xml_find_all(xmldoc,xmlpath)
  with(e,{
    eval(parse(text=paste("years=c(",xml_attr(detailed,"years"),")",sep="")))
    eval(parse(text=paste("age.from=c(",xml_attr(detailed,"age.group.lower.limits"),")",sep="")))
    combine_stocks_and_flows = xml_attr(xml_find_all(xmldoc,"//output/detailed.output"),"combine.stocks.and.flows")=="true"
    combine_stocks_and_flows = ifelse(is.na(combine_stocks_and_flows),F,combine_stocks_and_flows)
    suppress_zeros_flows = xml_attr(xml_find_all(xmldoc,"//output/detailed.output"),"suppress.zeros.in.flows")=="true"
    suppress_zeros_flows = ifelse(is.na(suppress_zeros_flows),F,suppress_zeros_flows)
    suppress_zeros_stocks = xml_attr(xml_find_all(xmldoc,"//output/detailed.output"),"suppress.zeros.in.stocks")=="true"
    suppress_zeros_stocks = ifelse(is.na(suppress_zeros_stocks),F,suppress_zeros_stocks)
    econ.output = xml_attr(xml_find_all(xmldoc,"//output/detailed.output"),"econ.output")=="true"
    econ.output = ifelse(is.na(econ.output),F,econ.output)
    options = output.options.from.xpath(xmldoc,"//output/detailed.output/flows")
  })
  e
}

#' a helper function to parse the <output><detailed.output><flows> elements
#' 
#' @param xmldoc the XML document
#' @param xpath by default "//output/detailed.output/flows"
#' 
#' @return an data frame with output options by dimension 
output.options.from.xpath=function(xmldoc,xpath){
  nodeset = xml_find_all(xmldoc,xpath)
  n=length(nodeset)
  if (n>0){
    Tm = rep(F,n) ; Tm[xml_attr(nodeset,"transmission")=="true"]=T
    Pr = rep(F,n) ; Pr[xml_attr(nodeset,"progression")=="true"]=T
    Ix = rep(F,n) ; Ix[xml_attr(nodeset,"incidence")=="true"]=T
    Tr = rep(F,n) ; Tr[xml_attr(nodeset,"treatment")=="true"]=T
    return(data.frame(dim=xml_attr(nodeset,"dim"), transmission=Tm, progression=Pr, incidence=Ix, treatment=Tr))
  }else{
    return(NULL)
  }
}

#' parses the incidence data specification for a dimension
#' 
#' @param xml   the XML document
#' @param xpath the path to the incidence specification in the XML e.g. <VXa><VXa.incidence>
#' 
#' @return a data frame with specification of incidence data ; columns file, times, values, proportions, denominator, once.per.year
inci.data.frame.from.xpath=function(xml,xpath){
  nodeset = xml_find_all(xml$doc,xpath)
  if (length(nodeset)>0){
    props=xml_attr(nodeset,"proportions")
    proportions = rep(F,length(props))
    proportions[props=="true"]=T
    oncepyr = xml_attr(nodeset,"once.per.year")
    once.per.year = rep(F,length(oncepyr))
    once.per.year[oncepyr=="true"]=T
    targetprev = xml_attr(nodeset,"target.prevalence")
    target.prevalence = rep(F,length(targetprev))
    target.prevalence[targetprev=="true"]=T    
    return(data.frame(file=xml_attr(nodeset,"file"),
                      times=xml_attr(nodeset,"times"),
                      values=xml_attr(nodeset,"values"),
                      proportions=proportions,
                      denominator=xml_attr(nodeset,"denominator"),
                      once.per.year=once.per.year,
                      target.prevalence=target.prevalence))
  }else{
    return(NULL)
  }
}
# not used
#parse.incidence.output.spec=function(xmldoc){
#  xmlpath = "//output/incidence.output"
#  if (length(xml_find_all(xmldoc,xmlpath))<=0){
#    return(NULL)
#  }  
#  incidence.output   = xml_find_all(xmldoc,xmlpath)
#  
#  include.output = xml_find_all(xmldoc,paste(xmlpath,"/include/output.dim",sep=""))
#  df = data.frame(dim          = xml_attr(include.output,"name"),
#                  transmission = xml_attr(include.output,"transmission"),
#                  progression  = xml_attr(include.output,"progression"),
#                  treatment    = xml_attr(include.output,"treatment"),
#                  from.data    = xml_attr(include.output,"incidence.from.data"),stringsAsFactors = F)
#  df[ is.na(df)]  = F
#  df[df=="true"]  = T
#  df[df=="false"] = F
#  df[df$dim!="TB",c("transmission","treatment")]=F
#  e = new.env()
#  with(e,{
#    eval(parse(text=paste("years=c(",xml_attr(incidence.output,"years"),")",sep="")))
#    output.dims = df
#  })
#  e
#}

