run_param_set <- function(model, cc, params, params_uid, vx_chars) {
  
  combined_ipj <- list()
  
  # set.paths = initialize the params
  model_paths <- model$set.paths(countrycode  = cc,
                                 xml          = vx_chars$xml,
                                 parameters   = vx_chars$input)
  
  if (grepl("baseline", vx_chars$runtype)){
    scen_baseline = NULL
  } else {
    scen_baseline = model$baseline_output 
  } 
  
  
  # Run the model with the parameter set
  output = model$run(model_paths, new.parameter.values = params, baseline = scen_baseline, output.flows = F)
  
  counts <- output$stocks
  counts <- counts[age_from == 0  & age_thru == 14, AgeGrp := "[0,14]"]
  counts <- counts[age_from == 15 & age_thru == 99, AgeGrp := "[15,99]"]
  counts <- counts[age_from == 0  & age_thru == 99, AgeGrp := "[0,99]"]
  counts <- counts[, !c("age_from", "age_thru")]
  
  #### Incidence
  inc <- counts[TB == "Dscount" & (HIV != "HIVdead"),]
  inc <- inc[, .(N_inc = sum(value)), by = .(Country = country, Year = year, AgeGrp, VXa, BMI = RISK)]
  inc <- inc[, Year := floor(Year)]
  
  #### Mortality
  mort <- counts[TB == "TBdead" & (HIV != "HIVdead"),]
  mort <- mort[, .(N_mort = sum(value)), by = .(Country = country, Year = year, AgeGrp, VXa, BMI = RISK)]
  mort <- mort[, Year := floor(Year)]
  
  #### Mortality - on treatment
  mort_tx <- counts[TB == "TTBdeadcount" & (HIV != "HIVdead"),]
  mort_tx <- mort_tx[, .(N_mort_tx = sum(value)), by = .(Country = country, Year = year, AgeGrp, VXa, BMI = RISK)]
  mort_tx <- mort_tx[, Year := floor(Year)]
  
  bmi <- counts[TB == "Un" | TB == "L0" | TB == "Ls" | TB == "Lf" | TB == "Dc" | TB == "Ds" | TB == "OT" | TB == "T" | TB == "R"]
  bmi <- bmi[, .(N_bmi = sum(value)), by = .(Country = country, Year = year, AgeGrp, VXa, BMI = RISK)]
  bmi <- bmi[, Year := floor(Year)]
  
  bmi_prev <- counts[TB == "Dc" | TB == "Ds"]
  bmi_prev <- bmi_prev[, .(N_bmi_prev = sum(value)), by = .(Country = country, Year = year, AgeGrp, VXa, BMI = RISK)]
  bmi_prev <- bmi_prev[, Year := floor(Year)]
  
  rm(counts)
  
  # Combine everything into one dataset
  n_epi <- bmi[inc, on = .(Country = Country, Year = Year, AgeGrp = AgeGrp, VXa = VXa, BMI = BMI)]
  n_epi <- n_epi[bmi_prev, on = .(Country = Country, Year = Year, AgeGrp = AgeGrp, VXa = VXa, BMI = BMI)]
  n_epi <- n_epi[mort, on = .(Country = Country, Year = Year, AgeGrp = AgeGrp, VXa = VXa, BMI = BMI)]
  n_epi <- n_epi[mort_tx, on = .(Country = Country, Year = Year, AgeGrp = AgeGrp, VXa = VXa, BMI = BMI)]
  
  rm(bmi)
  rm(bmi_prev)
  rm(inc)
  rm(mort)
  
  if (grepl("baseline", vx_chars$runtype)){
    model$baseline_output <- output
  }
  
  rm(output)
  
  # Add the vaccine characteristics
  combined_ipj[["n_epi"]] <- n_epi[, `:=`(uid     = params_uid,
                                          runtype = vx_chars$runtype)]
  combined_ipj
  
}

