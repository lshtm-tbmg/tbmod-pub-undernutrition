model_version_ = "3_3_2_E"
#' Returns the model version as a string
#' 
#' 
#' @return a string with the model version e.g. "3_3_2_E"
#' @export
modelversion = function(){
  return(model_version_)
}

schema_version_ = "TB-Vx-schema-V.xsd"
#' Returns the XML Schema version as a string
#' 
#' 
#' @return a string with XML Schema version e.g. "TB-Vx-schema-T.xsd"
#' 
#' @export 
schemaversion = function(){
  return(schema_version_)
}

#set.options=function(){
#  options("warnPartialMatchDollar"=T)
#}
#redirect.errors=function(paths){
#  sink(file(paths[["logfile"]],'wt'),append=T,type="message")
#}

#' Sets up the model log in the logs subfolder of the country folder as specified in the paths environment ; called from set.paths()
#'
#' @param paths   the paths environment from which the country folder will be derived
#' @param lglevel logging level
#' 
#' @return a logging function to be assigned to the global function modlog() that will be called for logging e.g. modlog(level="FATAL", msg="not good')
#' @export
setup.model.log=function(paths,lglevel="DEBUG"){
  assert_that(!is.na(paths[["countries"]]),msg="paths$countries should not be NULL")
  assert_that(!is.na(paths[["country.code"]]),msg="paths$country.code should not be NULL")
  RUNID_TS = Sys.getenv("RUNID_TS")
  if (RUNID_TS == "") {
    LTS <- format.Date(Sys.time(), format = "%Y-%m-%d-%H%M", tz = "UTC")
    RUNID_TS <- sprintf("%s_%s_%s", LTS, "LOCAL", paths$country.code)
  }
  logr <- create.logger(logfile = here(paths$log.dir, paste0(RUNID_TS, "_model.log")), level = lglevel)
  return(function(level = level, msg = NULL) { levellog(logr, level = level, message = msg)})
}

#' update.paths() creates an R environment with file paths to the XML file, parameter file, targets file etc
#'
#' For details see the set.paths() function that calls update.paths 
#'
#' 
#' @param paths the paths environment to be updated
#' @param mydir the ‘root’ directory
#' @param countries the name of the countries folder (which will be appended to mydir)
#' @param countrycode 3 letter ISO code of the country (e.g. ZAF)
#' @param xml the filename of the XML input file
#' @param parameters the filename of the parameters file (which may overwrite parameters in the xml inputfile)
#' @param targets the filename of the targets.csv file
#' @param baseline the filename (without extent) of the XML file of the baseline
#' @param VXa.incidence.files a vector of filenames (including path from the country folder) to replace the VXa incidence files specified in the XMLinput.xml file
#' @return the updated paths environment 
#' @export
update.paths = function(paths         = NULL,
                        mydir         = here(), 
                        countries     = "countries", 
                        countrycode   = NA, 
                        targets       = NA, 
                        xml           = NA, 
                        baseline      = NA, 
                        parameters    = NA, 
                        VXa.incidence.files=NA){
  assert_that(!is.null(paths),msg="paths should not be NULL ")
  assert_that(is.environment(paths),msg="paths should be an R environment")
  if (!(mydir)==here()){ paths$mydir = mydir}
  paths$country.code  = countrycode   
  paths$xml           = xml
  paths$parameters    = parameters
  paths$targets       = targets
  paths$baseline      = baseline
  paths$VXa.incidence.files = VXa.incidence.files 
  assert_that(!is.na(paths$country.code),msg="country code should not be NA")
  assert_that(!is.na(paths$xml),msg="xml should not be NA")
  paths$countries      = countries 
  paths$country.dir    = here(mydir,countries,countrycode) 
  paths$log.dir        = here(paths$country.dir,"logs")
  paths$data.dir       = here(paths$country.dir,"data")
  paths$params.dir     = here(paths$country.dir,"parameters")
  paths$output.dir     = here(paths$country.dir,"output")
  if (!is.na(paths$xml))       { paths$xml        = here(paths$params.dir,xml)}
  if (!is.na(paths$parameters)){ paths$parameters = here(paths$params.dir,parameters)}
  if (!is.na(paths$targets))   { paths$targets    = here(paths$params.dir,targets)}
  if (!is.na(paths$baseline))  { paths$baseline   = here(paths$params.dir,baseline)}
  if (!is.na(baseline)){
    popadj.fnames=dir(paths$output.dir,pattern=paste0("^\\[",countrycode,"\\].*\\[",baseline,"\\]\\[dfrPOPadj\\]\\.txt"))
    if (length(popadj.fnames)!=1){
      modlog(level='FATAL', msg= paste0("either none or multiple matches for baseline popadj output of ",baseline))
    }
    assert_that(length(popadj.fnames)==1,msg=paste0("either none or multiple matches for baseline popadj output of ",baseline))
    paths$popadj = here::here(paths$output.dir,popadj.fnames[1])
    bgxfr.fnames=dir(paths$output.dir,pattern=paste0("^\\[",countrycode,"\\].*\\[",baseline,"\\]\\[dfrBGx\\]\\.txt"))
    if (length(bgxfr.fnames)!=1){
      modlog(level='FATAL', msg= paste0("either none or multiple matches for baseline background death fractions output of ",baseline))
    }
    assert_that(length(bgxfr.fnames)==1,msg=paste0("either none or multiple matches for baseline background death fractions output of ",baseline))
    paths$bgxfr = here::here(paths$output.dir,bgxfr.fnames[1])
  }
  paths
}

#' set.paths() creates an R environment with file paths to the XML file, parameter file, targets file etc
#'
#' The only required parameters are countrycode and xml
#'
#' This function also initializes the model log
#' 
#' If the parameters argument is specified, the values of the parameters designated as constant in this file (usually input.csv) will overwrite the values in the XML.
#'
#' If the targets argument is specified, the output will contain a hits element with a comparison between targets and model values.
#
#' If the baseline argument is specified, the model run will use the adjustments that were made to the population size over time (and written to files) in the current simulation run to take into account reduced mortality in an intervention run with TB vaccines.
#' This correction requires the attribute value econ.output="true" in the XML inputfile (in <output><detailed.output>).
#'
#' The optional argument VXa.incidence.files is a vector of filenames that will replace the file names in <VXa><VXa.incidence><incidence.data>. Note: only the filenames will be replaced, not the additional attributes. Use with care.
#'
#' The directory structure of a specific country (ZAF as an example) and the required files are:
#' ZAF/data/demographics.csv 
#' ZAF/data/ZAF_deathrates.csv
#' ZAF/data/HIV-incidence.txt
#' ZAF/data/ART-incidence.txt
#' ZAF/data/ZAF_all_contacts_2020.txt
#' NOTE: the path and filename from the country folder depends on the specified path and filename in the XML input file e.g. <population file="data/demographics.csv"/>
#' ZAF/parameters/XMLinput.xml
#' ZAF/parameters/input.csv
#' ZAF/parameters/target.csv
#' NOTE: the filenames depend on the parameters passed to the intialize() function. The path is fixed.
#' ZAF/output	The folder where output will be written to
#' ZAF/logs       The folder where the model log will be written to
#' The initialize function also:
#' -	sets up the paths environment
#' -	sets up the global function modlog()
#' 
#' @param mydir the ‘root’ directory
#' @param countries the name of the countries folder (which will be appended to mydir)
#' @param countrycode 3 letter ISO code of the country (e.g. ZAF)
#' @param xml the filename of the XML input file
#' @param parameters the filename of the parameters file (which may overwrite parameters in the xml inputfile)
#' @param targets the filename of the targets.csv file
#' @param baseline the filename (without extent) of the XML file of the baseline
#' @param VXa.incidence.files a vector of filenames (including path from the country folder) to replace the VXa incidence files specified in the XMLinput.xml file
#' @return an environment with the paths to the various files and folders.
#' @export
set.paths = function(mydir         = here(), 
                     countries     = "countries", 
                     countrycode   = NA, 
                     xml           = NA, 
                     parameters    = NA, 
                     targets       = NA, 
                     baseline      = NA, 
                     VXa.incidence.files=NA,
                     lglevel       = "DEBUG",
                     taskenvvar    = "TaskID"){
  paths = new.env()
  paths$src = here("R")
  paths        = update.paths(paths,mydir,countries,countrycode,targets,xml,baseline,parameters,VXa.incidence.files)
  modlog       <<- setup.model.log(paths,lglevel)
  taskenvvar   <<- taskenvvar
  paths
}

#' Modifies the csv file with parameters that will replace parameters in the XML 
#' 
#' @param input.csv a data frame read from the file indicated by paths$parameters
#' @param new.parameter.values a named vector of parameter values to replace the non constant parameters in input.csv
#' @return the modified input.csv
#' @export
modify.input.csv = function(input.csv=NULL,new.parameter.values=NULL){
  assert_that(!is.null(input.csv),msg="input.csv argument should not be NULL in modify.input.csv()")  
  assert_that(!is.null(new.parameter.values),msg="new.parameter.values argument should not be NULL in modify.input.csv()")  
  constant   = constant.parameters(parameters.df=input.csv)
  err = any(constant$unique.name %in% names(new.parameter.values))
  assert_that(!err,msg="new.parameter.values should not contain parameter marked as constant")  
  if (err) {modlog(level="ERROR",msg="new.parameter.values should not contain parameter marked as constant")}
  sel = input.csv$unique.name %in% names(new.parameter.values)
  n   = nrow(input.csv)
  if (n>0){
   for (i in 1:nrow(input.csv[sel,])){
    input.csv[sel,][i,]$mean=new.parameter.values[input.csv[sel,][i,]$unique.name]
   }
  }
  return(input.csv)
}

#' Selects parameters from a data frame
#' 
#' @param parameters.df a data frame with parameters 
#' @param constant      a logical to indicate if constant or non constant parameters should be selected
#' @return a data frame with selected parameters
#' @export
selected.parameters = function(parameters.df=NULL, constant=NA){
  fncname = match.call()[[1]]
  assert_that(!is.null(parameters.df),msg=paste("parameters data frame should not be NULL in",fncname))
  assert_that(!is.na(constant),msg=paste("constant is a required argument [T/F] in",fncname))
  const            = grepl("^con.*", parameters.df$dist) | !parameters.df$choose
  if (constant){
    return(parameters.df[const,])
  }else{
    return(parameters.df[!const,])
  }  
}

#' Selects constant parameters from a data frame
#' 
#' @param parameters.df a data frame with parameters 
#' @return a data frame with constant parameters
#' @export
constant.parameters = function(parameters.df=NULL){
  return(selected.parameters(parameters.df,constant=T))
}

#' Selects non-constant parameters (to be fitted) from a data frame
#' 
#' @param parameters.df a data frame with parameters 
#' @return a data frame with non-constant parameters (to be fitted)
#' @export
fitted.parameters = function(parameters.df=NULL){
  return(selected.parameters(parameters.df,constant=F))
}

#' Samples uniform values from non-constant parameters from a data frame
#' 
#' @param selected.parameters a data frame with parameters 
#' @return a named vector with sampled parameter values (non-constant parameters only)
#' @export
sample.fitted.parameters = function(selected.parameters=NULL){
  fncname = match.call()[[1]]
  assert_that(!is.null(selected.parameters),msg=paste("selected.parameters should not be NULL in",fncname))
  constant            = grepl("^con.*", selected.parameters$dist) | !selected.parameters$choose
  varparams           = selected.parameters[!constant,]
  varparams$mean      = as.numeric(varparams$mean)
  varparams$min       = as.numeric(varparams$min)
  varparams$max       = as.numeric(varparams$max)
  newvalues           = runif(nrow(varparams),varparams$min,varparams$max)
  names(newvalues)    = varparams$unique.name
  newvalues
}

#' Reads a targets file, evaluates the output versus the targets and adds columns target hits to the targets file and returns a data frame
#' 
#' @param output model output
#' @param model.params environment with initialized model parameters
#' 
#' @return the model output with target hits added as a list element
#' @export
get.target.hits = function(output=NULL,model.params=NULL){
  assert_that(!is.null(model.params),msg="model.params should not be NULL when evaluating targets")
  if (!is.na(model.params$paths$targets)){
    targets     = read.targets(model.params$paths$targets)
    hits        = eval.output.vs.targets(output,targets,model.params)
    output$hits = cbind(country=model.params$run.params$countrycode,hits)
  }
  output
}

#' Writes the targets file with info on hit targets to a file
#' 
#' @param output model output
#' @param model.params environment with initialized model parameters
#' @param output.format file format
#' 
#' @return NULL
#' @export
write.targets = function(output=NULL,model.params=NULL, output.format='txt'){
  if(!is.null(output$hits)){
    outfname  = paste0(create.filename(model.params$paths$output.dir, model.params, runtype = "hits"),".txt")
    write.output(output$hits,outfname,output.format)
  }
}

#' Writes stocks and flows output
#' 
#' @param output model output
#' @param model.params environment with initialized model parameters
#' @param output.format file format
#' 
#' @return NULL
#' @export
write.stocks.and.flows = function(output=NULL,model.params=NULL, output.format='txt'){
  fncname = match.call()[[1]]
  assert_that(!is.null(model.params$paths),msg=paste("paths should not be NULL in",fncname))
  assert_that(!is.null(output),msg=paste("output should not be NULL in",fncname))
  filename = paste0(create.filename(model.params$paths$output.dir, model.params,runtype="stocks"),".txt")
  write.output(output$stocks,filename, output.format)
  #filename = paste0(create.filename(model.params$paths$output.dir, model.params,runtype="dead.999"),".txt")
  #write.output(output$dead.999,filename, output.format)
  #filename = paste0(create.filename(model.params$paths$output.dir, model.params,runtype="count.999"),".txt")
  #write.output(output$count.999,filename, output.format)
  #filename = paste0(create.filename(model.params$paths$output.dir, model.params,runtype="alive.999"),".txt")
  #write.output(output$alive.999,filename, output.format)
  #filename = paste0(create.filename(model.params$paths$output.dir, model.params,runtype="alive.500"),".txt")
  #write.output(output$alive.500,filename, output.format)
  if (!is.null(output$flows)){
    filename = paste0(create.filename(model.params$paths$output.dir, model.params,runtype="flows"),".txt")
    write.output(output$flows,filename, output.format)
  }
}

#' Writes econ output to files \code{###POP.txt}, \code{###dBGx.txt}, \code{###dHIVTBx.txt}, \code{###dPOPadj.txt}, \code{###dVx1pyr.txt}
#' 
#' @param output output from the run() function including econ output (population, dBGx, dHIVTBx, dPOPadj, vaccinated)     
#' @param model.params an environment with initialized model parameters
#' @param output.format the output format
#' 
#' @return NULL
#' @export
write.econ.output = function(output=NULL,model.params=NULL,output.format='txt',suffix=NA){
  
  fncname = match.call()[[1]]
  assert_that(!is.null(model.params),msg=paste("paths should not be NULL in",fncname))
  assert_that(!is.null(output),msg=paste("output should not be NULL in",fncname))
  
  if (!is.null(output$population)){
    filename = paste0(create.filename(model.params$paths$output.dir, model.params,runtype="POP"),".txt")
    write.output(output$population,filename,output.format)
    filename = paste0(create.filename(model.params$paths$output.dir, model.params,runtype="dBGx"),".txt")
    write.output(output$dBGx,filename,output.format)
    filename = paste0(create.filename(model.params$paths$output.dir, model.params,runtype="dTBHIVx"),".txt")
    write.output(output$dHIVTBx,filename,output.format)
    if (!model.params$intervention){
      filename = paste0(create.filename(model.params$paths$output.dir, model.params,runtype="dPOPadj"),".txt")
      write.output(output$dPOPadj,filename,output.format)
    }
    if (!is.null(output$vaccinated)){
      filename = paste0(create.filename(model.params$paths$output.dir, model.params,runtype="dVx1pyr"),".txt")
      write.output(output$vaccinated,filename,output.format)
    }
  }
  if (!is.na(suffix) & !is.null(output$population)){
    filename = paste0(create.filename(model.params$paths$output.dir, model.params,runtype="POP_NEW"),".txt")
    write.output(output$population.new,filename,output.format)
    filename = paste0(create.filename(model.params$paths$output.dir, model.params,runtype="dBGx_NEW"),".txt")
    write.output(output$dBGx.new,filename,output.format)
    filename = paste0(create.filename(model.params$paths$output.dir, model.params,runtype="dTBHIVx_NEW"),".txt")
    write.output(output$dHIVTBx.new,filename,output.format)
    if (!model.params$intervention){
      filename = paste0(create.filename(model.params$paths$output.dir, model.params,runtype="dPOPadj_NEW"),".txt")
      write.output(output$dPOPadj.new,filename,output.format)
    }
    if (!is.null(output$vaccinated)){
      filename = paste0(create.filename(model.params$paths$output.dir, model.params,runtype="dVx1pyr_NEW"),".txt")
      write.output(output$vaccinated.new,filename,output.format)
    }
  }
}

#' Writes the in memory XML doc that was updated with new parameter values to a file 
#' 
#' @param parameters environment with initialized (and updated) parameters
#' @param xml the new XML file name
#'  
#' @return NULL
#' @export
write.updated.xml = function(parameters=NULL, xml=NULL){
  assert_that(!is.null(parameters))
  assert_that(!is.null(xml))
  write_xml(parameters$xml$doc, file = here(parameters$paths$params.dir, xml), option="as_xml")
}

#' Merges stocks and flows output
#' 
#' @param output model output
#' 
#' @return data frame with merged stocks and flows model output
#' @export
merge.stocks.and.flows=function(output){
  widened.stocks = cbind(output$stocks[,1:7],dim=NA,subject=NA,flow=NA,output$stocks[,8:10])
  rbind(widened.stocks,output$flows)
}

#' Writes merged stocks and flows output
#' 
#' @param onerun merged stocks and flows model output
#' @param model.params environment with initialized model parameters
#' @param output.format file format
#' 
#' @return NULL
#' @export
write.merged.stocks.and.flows = function(onerun=NULL,model.params=NULL,output.format='txt'){
  fncname = match.call()[[1]]
  assert_that(!is.null(model.params$paths),msg=paste("paths should not be NULL in",fncname))
  assert_that(!is.null(onerun),msg=paste("onerun should not be NULL in",fncname))
  filename = paste0(create.filename(model.params$paths$output.dir, model.params,runtype="onerun"),".txt")
  write.output(onerun,filename,output.format)
}
