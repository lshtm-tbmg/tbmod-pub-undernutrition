#----------------------------------------------------------------------
# Generate Output for "Impacts of nutritional shocks on TB in India" paper
# Rebecca Clark
# Updated 6 February 2024
#----------------------------------------------------------------------

# 1. Load in the required packages
suppressPackageStartupMessages({
  rm(list=ls())
  model = new.env()
  library(here)
  library(data.table)
  library(renv)
  library(digest)
  library(log4r)
  library(fst)
  library(arrow)
  library(logger)
  
  source(here("R", "include-v11.R"), model)
  source(here("R", "TBVx-run-v1.R"), model)
  source(here("R", "workflow", "run_param_set.R"))
})


# 2. Set the country code, parameters, scenario characteristics 
cc <- "INDu"
grid_task_int <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

shock_chars_all  <- fread("./processing_files/INDu_shocks.csv")
scenarios_all    <- shock_chars_all[, .(scenario)]
scenarios_all    <- unique(scenarios_all)

parameters <- fread(paste0("./processing_files/param_sets/INDu_params.csv"))
parameters <- parameters[(10*(grid_task_int-1) + 1):(10*grid_task_int), ]
print(paste0("number of parameter sets to run = ", nrow(parameters)))


# 3. Set-up and run the model for each parameter set 

for (j in 1:nrow(parameters)) {
  
  print(paste0("parameter set = ", j))
  
  # select row j of parameters
  params <- parameters[j, ]
  
  # save the uid
  params_uid <- params[, uid]
  
  # get rid of everything except parameters
  params <- params[, !c("uid", "nhits")]
  params <- unlist(params)
  
  scenarios_j <- scenarios_all[sample(.N, 1)]
  shock_scenarios <- shock_chars_all[scenario == scenarios_j$scenario]
  print(vx_scenarios$runtype)
  
  cc_n_epi_param <- list()
  
  for (i in 1:nrow(vx_scenarios)){ # For each of the scenarios
    
    shock_chars <- shock_scenarios[i,]
    
    print(paste0("Running scenario number ", i, ": ", shock_chars$runtype))

    # run the model with the row of parameters
    shock_scen_output <- run_param_set(model, cc, params, params_uid, shock_chars)
    
    cc_n_epi_param[[i]] <- shock_scen_output[["n_epi"]]
  
    }

  write_parquet(rbindlist(cc_n_epi_param), here("epi_output_UKR", "n_epi", paste0(params_uid, ".parquet")))
  rm(cc_n_epi_param)
  
  print(paste0("End time for parameter set ", j, " = ", Sys.time()))

  }


# ----end

