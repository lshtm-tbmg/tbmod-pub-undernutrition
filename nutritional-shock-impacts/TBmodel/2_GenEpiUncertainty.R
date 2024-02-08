#-------------------
# Generate uncertainty bounds for output generated from 1_GenEpiOutput.R
# Rebecca Clark
# Last updated 8 February 2024
#-------------------

# 0. Set-up
suppressPackageStartupMessages({
  rm(list=ls())
  library(here)
  library(data.table)
  library(renv)
  library(arrow)
  library(magrittr)
})

# Load in data generated from GenEpiOutput.R
n_epi <- open_dataset(sources = "./epi_output_UKR/n_epi/") %>% dplyr::collect()
n_epi <- setDT(n_epi)

n_epi <- n_epi[grepl("baseline", n_epi$runtype), Scenario := "No Shocks"]
n_epi <- n_epi[grepl("energy_export", n_epi$runtype), Scenario := "Export and Energy Shocks"]
n_epi <- n_epi[is.na(Scenario) & grepl("export_", n_epi$runtype), Scenario := "Export Restriction Shock"]
n_epi <- n_epi[is.na(Scenario) & grepl("energy", n_epi$runtype), Scenario := "Energy Price Shock"]


#----------------------------------------------------------------------
# 1. Calculating incidence/mortality rates and number of cases/deaths/deaths on treatment
#----------------------------------------------------------------------

total <- n_epi[,`:=`(N_inc = sum(N_inc),
                     N_mort = sum(N_mort),
                     N_mort_tx = sum(N_mort_tx),
                     N_total = sum(N_bmi)),
               by = .(Year, AgeGrp, Scenario, BMI, uid)]

total <- total[, .(Year, AgeGrp, BMI, Scenario, UID = uid, N_inc, N_mort, N_mort_tx, N_total)]
total <- unique(total)

total_raw_rate <- total[, `:=`(inc_rate = (N_inc/N_total)*(100000),
                               mort_rate = (N_mort/N_total)*(100000),
                               mort_tx_rate = (N_mort_tx/N_total)*(100000)),
                        by = .(Year, AgeGrp, Scenario, UID, BMI)]

# sum = number of cases/deaths/deaths on tx per year by scenario
total_raw_sum <- total[, `:=`(sum_inc = sum(N_inc),
                              sum_mort = sum(N_mort),
                              sum_mort_tx = sum(N_mort_tx)),
                       by = .(Year, AgeGrp, Scenario, UID, BMI)]

total_raw <- merge(total_raw_sum, total_raw_rate)

total_raw <- total_raw[, .(UID, Year, AgeGrp, Scenario, BMI,
                           N_inc, N_mort, N_mort_tx, N_total, inc_rate, mort_rate,
                           mort_tx_rate, sum_inc, sum_mort, sum_mort_tx)]

total_raw_long <- melt(total_raw, id.vars = c("Year", "AgeGrp", "Scenario", "UID", "BMI"),
                       measure.vars = c("inc_rate", "mort_rate", "mort_tx_rate",
                                        "sum_inc", "sum_mort", "sum_mort_tx"),
                       variable.name = "Indicator", value.name = "Value")


total_raw_long <- setDT(total_raw_long)
total_raw_long <- total_raw_long[, .(medval = median(Value),
                                     lowval = quantile(Value, 0.025),
                                     highval = quantile(Value, 0.975)),
                                 by = .(Year, AgeGrp, Scenario, BMI, Indicator)]

total_raw_long <- total_raw_long[, combined := paste0(round(medval, 1),
                                                      " (", round(lowval, 1),
                                                      ", ", round(highval, 1), ")")]


fwrite(total_raw_long, "./epi_output_UKR/raw_output_undernutrition_UKR_BMI.csv")

# Cumulative number of cases/deaths from a starting year
total_raw_cumulative <- total_raw_sum[Year >= 2021,]

total_raw_cumulative <- total_raw_cumulative[, `:=`(sum_inc = sum(N_inc),
                                                    sum_mort = sum(N_mort),
                                                    sum_mort_tx = sum(N_mort_tx)),
                                             by = .(Year, AgeGrp, Scenario, UID)]

total_raw_cumulative <- total_raw_cumulative[, .(Year, AgeGrp, Scenario, UID, sum_inc, sum_mort, sum_mort_tx)]

total_raw_cumulative <- unique(total_raw_cumulative)

total_raw_cumulative <- total_raw_cumulative[, `:=`(sum_inc = cumsum(sum_inc),
                                                    sum_mort = cumsum(sum_mort),
                                                    sum_mort_tx = cumsum(sum_mort_tx)),
                                             by = .(AgeGrp, Scenario, UID)]

total_raw_cumulative <- total_raw_cumulative[, .(Year, AgeGrp, Scenario, UID, sum_inc, sum_mort, sum_mort_tx)]
total_raw_cumulative <- unique(total_raw_cumulative)


#----------------------------------------------------------------------
# 2. Calculating incidence and mortality rate reductions/increases (compared to No Shocks)
#----------------------------------------------------------------------

rate_red_long <- melt(total_raw_rate, measure.vars = c("inc_rate", "mort_rate"), 
                      id.vars = c("Year", "AgeGrp", "Scenario", "UID", "BMI"),
                      value.name = "Value", variable.name = "Indicator")

rate_red_wide <- dcast(rate_red_long, Year + UID + BMI + AgeGrp + Indicator ~ Scenario,
                       value.var = "Value")

shock_scenarios <- unique(total_raw_rate$Scenario)
shock_scenarios <- shock_scenarios[shock_scenarios != "No Shocks"]

rate_red_wide <- setDT(rate_red_wide)

for (scen in shock_scenarios) {
  set(x = rate_red_wide, j = paste0(scen, "_diff"), value = rate_red_wide[["No Shocks"]] - rate_red_wide[[paste0(scen)]])
  set(x = rate_red_wide, j = paste0("PER_", scen), value = rate_red_wide[[paste0(scen, "_diff")]] / rate_red_wide[["No Shocks"]])
}

rate_reductions <- melt(data = rate_red_wide, 
                        measure.vars = c("PER_Energy Price Shock", "PER_Export Restriction Shock", "PER_Export and Energy Shocks"), 
                        id.vars = c("UID", "AgeGrp", "Year", "BMI", "Indicator"),
                        value.name = "Value", variable.name = "Scenario")

rate_reductions <- setDT(rate_reductions)
rate_reductions <- rate_reductions[Indicator == "inc_rate", Indicator := "inc_RR"]
rate_reductions <- rate_reductions[Indicator == "mort_rate", Indicator := "mort_RR"]

rate_reductions <- rate_reductions[, .(medval = median(Value),
                                       lowval = quantile(Value, 0.025),
                                       highval = quantile(Value, 0.975)),
                                   by = .(Year, Scenario, AgeGrp, BMI, Indicator)]

rate_reductions <- rate_reductions[, combined := paste0(round(medval*100, 1),
                                                        "% (", round(lowval*100, 1),
                                                        ", ", round(highval*100, 1), ")")]

rate_reductions$Scenario <- gsub("PER_", "", as.character(rate_reductions$Scenario))

fwrite(rate_reductions, "./epi_output_UKR/rate_differences_UKR_BMI.csv")



#-------------------------------------------------------------------------------------------
# 3.  Calculating cumulative tx, cases, and deaths averted/added by specific years
#-------------------------------------------------------------------------------------------

sum_epi_long <- melt(total_raw_cumulative, measure.vars = c("sum_inc", "sum_mort"), 
                     id.vars = c("Year", "AgeGrp", "Scenario", "UID"),
                     value.name = "Value", variable.name = "Indicator")

sum_epi_wide <- dcast(sum_epi_long, Year + UID + AgeGrp + Indicator ~ Scenario,
                      value.var = "Value")

shock_scenarios <- unique(sum_epi_long$Scenario)
shock_scenarios <- shock_scenarios[shock_scenarios  != "No Shocks"]

for (scen in shock_scenarios) {
  set(x = sum_epi_wide, j = paste0("Diff_", scen), value = sum_epi_wide[["No Shocks"]] - sum_epi_wide[[paste0(scen)]])
}

sum_averted <- melt(data = sum_epi_wide, measure.vars = patterns("^Diff.*"), 
                    id.vars = c("Year", "UID", "AgeGrp", "Indicator"),
                    value.name = "Value", variable.name = "Scenario")

sum_averted$Runtype <- gsub("Diff_", "", as.character(sum_averted$Runtype))

sum_averted <- sum_averted[Indicator == "sum_inc",  Indicator := "inc_avert"]
sum_averted <- sum_averted[Indicator == "sum_mort", Indicator := "mort_avert"]
sum_averted <- sum_averted[, .(medval = median(Value),
                               lowval = quantile(Value, 0.025),
                               highval = quantile(Value, 0.975)),
                           by = .(Year, Scenario, Indicator, AgeGrp)]

sum_averted <- sum_averted[, combined := paste0(round(medval, 1),
                                                " (", round(lowval, 1),
                                                ", ", round(highval, 1), ")")]






# ---- end


