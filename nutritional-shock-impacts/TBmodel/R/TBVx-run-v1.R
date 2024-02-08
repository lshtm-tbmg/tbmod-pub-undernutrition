#'The purpose of run() is to read the XML input file, update parameter values (if applicable), intialize the model parameters (i.e. setting up the matrices), do one model run, and return the output.
#' 
#' If the argument new.parameter.values is not supplied while the paths$parameters is not NA, all parameters from the parameters file (input.csv) will be used to overwrite the parameter values in the XML.
#' 
#' @param paths=NULL initialized paths
#' @param new.parameter.values=NULL a named vector of parameter values that will overwrite (part of) the non-constant parameters in the parameters file (usually input.csv)
#' @param write.to.file=F a logical that determines if output will be written to files
#' @param write.xml=NULL the name of a file to write the modified XML to
#' @param output.flows=T a logical that determines flows will be calculated and included in the output
#' @param combine.stocks.and.flows=F if T stocks and flows will be written to a single file 
#' @param baseline=NULL if not NULL, a reference to an output structure produced by a baseline run
#' @param output.format='txt' either ‘txt’ for tab-delimited text, ‘fst’ for fst, or ‘parquet’ for parquet
#' @param sample.parameters=F if T, values of the non-constant parameters in input,csv will be sampled from the distribution specified in the ‘distr’ column of input.csv 
#' 
#' 
#' @return a list with the following elements: stocks, flows, hits (only if targets is specified) and additional elements if econ.output=”true” in the XML
#' 
#' @export
run=function(paths=NULL,new.parameter.values=NULL,write.to.file=F,write.xml=NULL, output.flows=T,
             combine.stocks.and.flows=F, baseline=NULL, output.format='txt', sample.parameters = F){

  input.csv = NULL
  
  if (is.na(paths$parameters)){
    assert_that(is.null(new.parameter.values),msg="path to input.csv should not be NA when replacing parameter values with new.parameter.values")
  }else{
    input.csv = read.csv(paths$parameters, stringsAsFactors = F, header=T, fileEncoding = 'UTF-8-BOM')
    if (!is.null(new.parameter.values)){
      if(any(class(new.parameter.values)=="data.table")){
        new.parameter.values = as.data.frame(new.parameter.values)
      }
    }else if (sample.parameters){
      new.parameter.values = sample.fitted.parameters(selected.parameters = input.csv)
    }  
    if (!is.null(new.parameter.values)){
      input.csv  = modify.input.csv(input.csv,new.parameter.values)
    }
  }
  
  xmlparams  = read.model.parameters(paths)
  
  if (!is.null(baseline)){
    set.baseline.in.model.parameters(xmlparams,output=baseline)
  }

  if (!is.null(input.csv)){
    xmlparams$xml$doc = update.parameters(xmlparams$xml$doc,input.csv)
  }
  
  initialized.params = initialize.model.parameters(xmlparams)
  
  if (!is.null(write.xml)){
    write.updated.xml(initialized.params, xml=write.xml)
  }
  output = run.model(initialized.params, output.flows)
  gc()
  if (!is.null(output$stocks)){
   output  = get.target.hits(output, initialized.params)
   gc()
   if (write.to.file){
    if (combine.stocks.and.flows){
      onerun = merge.stocks.and.flows(output)
      write.merged.stocks.and.flows(onerun,initialized.params, output.format)
    }else{
      write.stocks.and.flows(output,initialized.params, output.format) 
    }
    write.targets(output, initialized.params, output.format) 
    write.econ.output(output, initialized.params, output.format) 
   }        
  }
  output
}
