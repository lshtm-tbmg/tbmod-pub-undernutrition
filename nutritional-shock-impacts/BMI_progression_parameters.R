## Finn McQuaid 15/12/2022
## Estimation of BMI-weighted progression and reactivation parameters
## Based on the assumption that infection rate is the same
## Using Incidence, BMI distribution, and the log-linear relationship between the two

## 1. Background and required packages =========================================
## We set the working file directory and load all required R packages  
setwd("C:/Users/eidefmcq/Documents/Simulations/Undernutrition")
rm(list=ls())
library("ggplot2")
library("data.table")
library("tidyverse")
library("rriskDistributions")
library("survey")
library("ipumsr")
library("rdhs")

## 2. Loading data =============================================================
## We load the required data
## TB incidence estimates from https://www.who.int/teams/global-tuberculosis-programme/data
## We use overall population estimates, not disaggregating by age or sex
inc_est   <- fread('TB_burden_countries_2022-12-02.csv')
inc_est   <- inc_est[country=="India"]

## BMI and self-reported TB from https://www.idhsdata.org/idhs/index.shtml (HWMBMI, HWFBMI and TBSUFFER)
## Excluding survey from 1998-1999 as men were unavailable
## Excluding 2019-2021 as TB incidence estimates were affected by COVID-19
## BMI distribution for younger adults (not older adults or children)
ddi       <- read_ipums_ddi("idhs_00004.xml")
BMI_data  <- read_ipums_micro(ddi)
BMI_data <- filter(BMI_data, BMI_data$YEAR==2015)
## Merging male and female results into one column for preparation
## Divide by 100 (BMI is recorded x100 to remove decimal places)
BMI_data  <- transform(BMI_data, bmi=pmin(HWMBMI, HWFBMI)/100)
## Remove rows with flagged or missing values
BMI_data <- filter(BMI_data, BMI_data$bmi<9997)
## Adding a categorical variable to indicate BMI category (<17.0, 17.0-18.5, 18.5-25, >25)
BMI_data$bmi_cat <- as.factor(ifelse(BMI_data$bmi<17.0, 'mod',
                                     ifelse(BMI_data$bmi<18.5, 'mild',
                                            ifelse(BMI_data$bmi<25, 'normal','over'))))
## Use survey to weight population-level results http://asdfree.com/demographic-and-health-surveys-dhs.html
dhs_design <- svydesign(~IDHSPSU, strata=~IDHSSTRATA, data=BMI_data, weights=~HHWEIGHT)
dhs_design <- update(dhs_design, one=1, bmi=bmi, tb=TBSUFFER, bmi_cat=bmi_cat, HWFBMI=HWFBMI, HWMBMI=HWMBMI)
## Removing strata have only 1 value
options(survey.adjust.domain.lonely=TRUE)
options(survey.lonely.psu="adjust")
## Separating out male and female for calculations
sub_dhs_design_M <- subset(dhs_design, HWMBMI<9997)
sub_dhs_design_F <- subset(dhs_design, HWFBMI<9997)
## Average BMI of those in each BMI category
BMI_avg_M    <- svyby(~bmi, ~bmi_cat, sub_dhs_design_M, svymean)
BMI_avg_F    <- svyby(~bmi, ~bmi_cat, sub_dhs_design_F, svymean)
## Proportion of population in each BMI category
BMI_pop_M    <- svyby(~one, ~bmi_cat, sub_dhs_design_M, svytotal)
BMI_pop_F    <- svyby(~one, ~bmi_cat, sub_dhs_design_F, svytotal)
## For the total population use F:M sex ratio from https://www.statcompiler.com/en/
FM_pop       <- c(FM_ratio=1, FM_ratio_lo=0.994, FM_ratio_hi=1.007)
## Get log normal distributions for the above from mean & CI for sampling later
modthin_pop_F_lnorm  <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(BMI_pop_F)[2,1],BMI_pop_F[2,2],confint(BMI_pop_F)[2,2]),show.output=FALSE,plot=FALSE)))
modthin_pop_M_lnorm  <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(BMI_pop_M)[2,1],BMI_pop_M[2,2],confint(BMI_pop_M)[2,2]),show.output=FALSE,plot=FALSE)))
mildthin_pop_F_lnorm <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(BMI_pop_F)[1,1],BMI_pop_F[1,2],confint(BMI_pop_F)[1,2]),show.output=FALSE,plot=FALSE)))
mildthin_pop_M_lnorm <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(BMI_pop_M)[1,1],BMI_pop_M[1,2],confint(BMI_pop_M)[1,2]),show.output=FALSE,plot=FALSE)))
normal_pop_F_lnorm   <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(BMI_pop_F)[3,1],BMI_pop_F[3,2],confint(BMI_pop_F)[3,2]),show.output=FALSE,plot=FALSE)))
normal_pop_M_lnorm   <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(BMI_pop_M)[3,1],BMI_pop_M[3,2],confint(BMI_pop_M)[3,2]),show.output=FALSE,plot=FALSE)))
over_pop_F_lnorm     <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(BMI_pop_F)[4,1],BMI_pop_F[4,2],confint(BMI_pop_F)[4,2]),show.output=FALSE,plot=FALSE)))
over_pop_M_lnorm     <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(BMI_pop_M)[4,1],BMI_pop_M[4,2],confint(BMI_pop_M)[4,2]),show.output=FALSE,plot=FALSE)))
modthin_avg_F_lnorm  <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(BMI_avg_F)[2,1],BMI_avg_F[2,2],confint(BMI_avg_F)[2,2]),show.output=FALSE,plot=FALSE)))
modthin_avg_M_lnorm  <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(BMI_avg_M)[2,1],BMI_avg_M[2,2],confint(BMI_avg_M)[2,2]),show.output=FALSE,plot=FALSE)))
mildthin_avg_F_lnorm <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(BMI_avg_F)[1,1],BMI_avg_F[1,2],confint(BMI_avg_F)[1,2]),show.output=FALSE,plot=FALSE)))
mildthin_avg_M_lnorm <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(BMI_avg_M)[1,1],BMI_avg_M[1,2],confint(BMI_avg_M)[1,2]),show.output=FALSE,plot=FALSE)))
normal_avg_F_lnorm   <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(BMI_avg_F)[3,1],BMI_avg_F[3,2],confint(BMI_avg_F)[3,2]),show.output=FALSE,plot=FALSE)))
normal_avg_M_lnorm   <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(BMI_avg_M)[3,1],BMI_avg_M[3,2],confint(BMI_avg_M)[3,2]),show.output=FALSE,plot=FALSE)))
over_avg_F_lnorm     <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(BMI_avg_F)[4,1],BMI_avg_F[4,2],confint(BMI_avg_F)[4,2]),show.output=FALSE,plot=FALSE)))
over_avg_M_lnorm     <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(BMI_avg_M)[4,1],BMI_avg_M[4,2],confint(BMI_avg_M)[4,2]),show.output=FALSE,plot=FALSE)))

## Estimating self-reported TB from DHS data by BMI category, subset with data only
sub_dhs_design_TB <- subset(dhs_design, TBSUFFER<2)
## Separating out male and female for calculations
sub_dhs_design_M_TB <- subset(sub_dhs_design_TB, HWMBMI<9997)
sub_dhs_design_F_TB <- subset(sub_dhs_design_TB, HWFBMI<9997)
## TB incidence per 100,000 for men and women
TB_M    <- svyby(~tb, ~bmi_cat, sub_dhs_design_M_TB, svymean)
TB_F    <- svyby(~tb, ~bmi_cat, sub_dhs_design_F_TB, svymean)
TB_F_modthin_lnorm  <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(TB_F)[2,1],TB_F[2,2],confint(TB_F)[2,2]),show.output=FALSE,plot=FALSE)))
TB_M_modthin_lnorm  <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(TB_M)[2,1],TB_M[2,2],confint(TB_M)[2,2]),show.output=FALSE,plot=FALSE)))
TB_F_mildthin_lnorm <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(TB_F)[1,1],TB_F[1,2],confint(TB_F)[1,2]),show.output=FALSE,plot=FALSE)))
TB_M_mildthin_lnorm <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(TB_M)[1,1],TB_M[1,2],confint(TB_M)[1,2]),show.output=FALSE,plot=FALSE)))
TB_F_normal_lnorm     <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(TB_F)[3,1],TB_F[3,2],confint(TB_F)[3,2]),show.output=FALSE,plot=FALSE)))
TB_M_normal_lnorm     <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(TB_M)[3,1],TB_M[3,2],confint(TB_M)[3,2]),show.output=FALSE,plot=FALSE)))
TB_F_over_lnorm     <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(TB_F)[4,1],TB_F[4,2],confint(TB_F)[4,2]),show.output=FALSE,plot=FALSE)))
TB_M_over_lnorm     <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(confint(TB_M)[4,1],TB_M[4,2],confint(TB_M)[4,2]),show.output=FALSE,plot=FALSE)))

## BMI/incidence relationship from https://pubmed.ncbi.nlm.nih.gov/19820104/
## Incidence for Cegielski in Lonroth is incorrect, use instead ttps://academic.oup.com/aje/article/176/5/409/100183?login=false 
## For Cegielski BMI average assume a uniform distribution where possible, and for <18 and >30 use the original values from Lonnroth
## Note lack of data for low BMI, such that the RR we find for underweight compartment is likely a lower bound
palmer    <- data.frame(author=rep("palmer",5),    bmi=c(18.5,19.8,21,22.6,25),   inc=c(75.1,67.3,35.2,29.7,18.9),      inc_lo=c(43.7,42.8,19.8,18.3,7.7),       inc_hi=c(106.5,91.8,50.7,41.1,30.0))
edwards   <- data.frame(author=rep("edwards",3),   bmi=c(18.5,21.8,25),           inc=c(79.1,46.5,22.5),                inc_lo=c(64.7,40.4,15.8),                inc_hi=c(93.6,52.6,29.1))
hemila    <- data.frame(author=rep("hemila",3),    bmi=c(21,25,29),               inc=c(116.0,55.2,32.7),               inc_lo=c(79,38.1,19.0),                  inc_hi=c(153,72.3,46.3))
tverdal   <- data.frame(author=rep("tverdal",7),   bmi=c(18.5,22,24,26,28,30,32), inc=c(16.7,11.2,8.3,6.2,4.1,2.8,3.1), inc_lo=c(15.3,10.2,7.4,5.4,3.3,1.8,2.1), inc_hi=c(18.0,12.1,9.1,7.0,4.9,3.7,4.0))
leung     <- data.frame(author=rep("leung",5),     bmi=c(17,20.75,24,27.5,31.5),  inc=c(599,291,218,163,102),           inc_lo=c(425,251,163,122,47),            inc_hi=c(847,335,294,218,218))
cegielski <- data.frame(author=rep("cegielski",4), bmi=c(17.5,21.75,27.5,34.2),   inc=c(260.2,24.7,8.9,5.1),            inc_lo=c(98.6,13.0,2.2,0.01),            inc_hi=c(421.8,36.3,15.6,10.5))
bmi_inc_rel <- rbind(palmer, edwards, hemila, tverdal, leung, cegielski)

## 3. Sampling data to use =====================================================
## We sample from the various distributions above to fit a linear model to the data
sample_sz <- 10000
rr_modthin  <- rep(1,sample_sz)
rr_mildthin <- rep(1,sample_sz)
rr_over     <- rep(1,sample_sz)
tb_modthin  <- rep(1,sample_sz)
tb_mildthin <- rep(1,sample_sz)
tb_normal <- rep(1,sample_sz)
tb_over     <- rep(1,sample_sz)
for (y in 1:sample_sz){
  ## Use either the 2005/06 or 2015/2016 data, here we've elected to use 2015
  #year_sample      <- sample(c(2005,2015), size = 1)
  year_sample       <- 2015
  ## Sample incidence from WHO estimates and 95% CI (which are the 2.5, 50 and 97.5 percentiles), log-normally distributed as it can't go negative
  inc_est_year     <- c(inc_est[year==year_sample,e_inc_100k_lo],inc_est[year==year_sample,e_inc_100k],inc_est[year==year_sample,e_inc_100k_hi])
  inc_est_lnorm    <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=inc_est_year,show.output=FALSE,plot=FALSE)))
  inc_est_sample   <- rlnorm(1, meanlog = inc_est_lnorm[1], sdlog = inc_est_lnorm[2])
  ## Sample for DHS for a sampled FM ratio
  #DHS_sample          <- DHS_bmi[DHS_bmi$year==year_sample,]
  female_modthin_pop  <- rlnorm(1, meanlog =  modthin_pop_F_lnorm[1], sdlog =  modthin_pop_F_lnorm[2])
  female_mildthin_pop <- rlnorm(1, meanlog=mildthin_pop_F_lnorm[1], sdlog=mildthin_pop_F_lnorm[2])
  female_normal_pop   <- rlnorm(1, meanlog=normal_pop_F_lnorm[1], sdlog=normal_pop_F_lnorm[2])
  female_over_pop     <- rlnorm(1, meanlog=over_pop_F_lnorm[1], sdlog=over_pop_F_lnorm[2])
  total_pop_F         <- female_modthin_pop + female_mildthin_pop + female_normal_pop  + female_over_pop
  male_modthin_pop    <- rlnorm(1, meanlog =  modthin_pop_M_lnorm[1], sdlog =  modthin_pop_M_lnorm[2])
  male_mildthin_pop   <- rlnorm(1, meanlog=mildthin_pop_M_lnorm[1], sdlog=mildthin_pop_M_lnorm[2])
  male_normal_pop     <- rlnorm(1, meanlog=normal_pop_M_lnorm[1], sdlog=normal_pop_M_lnorm[2])
  male_over_pop       <- rlnorm(1, meanlog=over_pop_M_lnorm[1], sdlog=over_pop_M_lnorm[2])
  total_pop_M         <- male_modthin_pop + male_mildthin_pop + male_normal_pop  + male_over_pop
  DHS_sample <- data.frame(female_modthin_prop  = female_modthin_pop/total_pop_F,
                         female_mildthin_prop   = female_mildthin_pop/total_pop_F,
                         female_normal_prop     =  female_normal_pop/total_pop_F,
                         female_over_prop       = female_over_pop/total_pop_F,
                         male_modthin_prop      = male_modthin_pop/total_pop_M,
                         male_mildthin_prop     = male_mildthin_pop/total_pop_M,
                         male_normal_prop       =  male_normal_pop/total_pop_M,
                         male_over_prop         = male_over_pop/total_pop_M,
                         female_modthin_avg  = rlnorm(1, meanlog =  modthin_avg_F_lnorm[1], sdlog =  modthin_avg_F_lnorm[2]),
                         female_mildthin_avg = rlnorm(1, meanlog =  mildthin_avg_F_lnorm[1], sdlog =  mildthin_avg_F_lnorm[2]),
                         female_normal_avg   = rlnorm(1, meanlog =  normal_avg_F_lnorm[1], sdlog =  normal_avg_F_lnorm[2]),
                         female_over_avg     = rlnorm(1, meanlog =  over_avg_F_lnorm[1], sdlog =  over_avg_F_lnorm[2]),
                         male_modthin_avg    = rlnorm(1, meanlog =  modthin_avg_M_lnorm[1], sdlog =  modthin_avg_M_lnorm[2]),
                         male_mildthin_avg   = rlnorm(1, meanlog =  mildthin_avg_M_lnorm[1], sdlog =  mildthin_avg_M_lnorm[2]),
                         male_normal_avg     = rlnorm(1, meanlog =  normal_avg_M_lnorm[1], sdlog =  normal_avg_M_lnorm[2]),
                         male_over_avg       = rlnorm(1, meanlog =  over_avg_M_lnorm[1], sdlog =  over_avg_M_lnorm[2]),
                         female_modthin_tb  = rlnorm(1, meanlog =  TB_F_modthin_lnorm[1], sdlog =  TB_F_modthin_lnorm[2]),
                         female_mildthin_tb = rlnorm(1, meanlog =  TB_F_mildthin_lnorm[1], sdlog =  TB_F_mildthin_lnorm[2]),
                         female_normal_tb   = rlnorm(1, meanlog =  TB_F_normal_lnorm[1], sdlog =  TB_F_normal_lnorm[2]),
                         female_over_tb     = rlnorm(1, meanlog =  TB_F_over_lnorm[1], sdlog =  TB_F_over_lnorm[2]),
                         male_modthin_tb    = rlnorm(1, meanlog =  TB_M_modthin_lnorm[1], sdlog =  TB_M_modthin_lnorm[2]),
                         male_mildthin_tb   = rlnorm(1, meanlog =  TB_M_mildthin_lnorm[1], sdlog =  TB_M_mildthin_lnorm[2]),
                         male_normal_tb     = rlnorm(1, meanlog =  TB_M_normal_lnorm[1], sdlog =  TB_M_normal_lnorm[2]),
                         male_over_tb       = rlnorm(1, meanlog =  TB_M_over_lnorm[1], sdlog =  TB_M_over_lnorm[2]))   
  FM_ratio_lnorm    <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(FM_pop[['FM_ratio_lo']],FM_pop[['FM_ratio']],FM_pop[['FM_ratio_hi']]),show.output=FALSE,plot=FALSE)))
  FM_ratio_sample   <- rlnorm(1, meanlog =  FM_ratio_lnorm[1], sdlog =  FM_ratio_lnorm[2])
  ## Calculate total proportion in different BMI categories 
  DHS_sample$modthin_prop  <- (DHS_sample[['male_modthin_prop']] + FM_ratio_sample*DHS_sample[['female_modthin_prop']]) /2
  DHS_sample$mildthin_prop <- (DHS_sample[['male_mildthin_prop']] + FM_ratio_sample*DHS_sample[['female_mildthin_prop']]) /2
  DHS_sample$normal_prop   <- (DHS_sample[['male_normal_prop']] + FM_ratio_sample*DHS_sample[['female_normal_prop']]) /2
  DHS_sample$over_prop     <- (DHS_sample[['male_over_prop']] + FM_ratio_sample*DHS_sample[['female_over_prop']]) /2
  ## Estimate BMI averages, weighting by FM ratio plus FM ratio for each separate BMI category
  DHS_sample$modthin_avg  <- (DHS_sample[['male_modthin_avg']] + FM_ratio_sample*DHS_sample[['female_modthin_avg']]*DHS_sample[['female_modthin_prop']]/DHS_sample[['male_modthin_prop']]) /2
  DHS_sample$mildthin_avg <- (DHS_sample[['male_mildthin_avg']] + FM_ratio_sample*DHS_sample[['female_mildthin_avg']]*DHS_sample[['female_mildthin_prop']]/DHS_sample[['male_mildthin_prop']]) /2
  DHS_sample$normal_avg   <- (DHS_sample[['male_normal_avg']] + FM_ratio_sample*DHS_sample[['female_normal_avg']]*DHS_sample[['female_normal_prop']]/DHS_sample[['male_normal_prop']]) /2
  DHS_sample$over_avg     <- (DHS_sample[['male_over_avg']] + FM_ratio_sample*DHS_sample[['female_over_avg']]*DHS_sample[['female_over_prop']]/DHS_sample[['male_over_prop']]) /2
  ## Estimate TB incidence averages, weighting by FM ratio plus FM ratio for each separate BMI category
  tb_modthin[y]  <- (DHS_sample[['male_modthin_tb']] + FM_ratio_sample*DHS_sample[['female_modthin_tb']]*DHS_sample[['female_modthin_prop']]/DHS_sample[['male_modthin_prop']]) /2
  tb_mildthin[y] <- (DHS_sample[['male_mildthin_tb']] + FM_ratio_sample*DHS_sample[['female_mildthin_tb']]*DHS_sample[['female_mildthin_prop']]/DHS_sample[['male_mildthin_prop']]) /2
  tb_normal[y]   <- (DHS_sample[['male_normal_tb']] + FM_ratio_sample*DHS_sample[['female_normal_tb']]*DHS_sample[['female_normal_prop']]/DHS_sample[['male_normal_prop']]) /2
  tb_over[y]     <- (DHS_sample[['male_over_tb']] + FM_ratio_sample*DHS_sample[['female_over_tb']]*DHS_sample[['female_over_prop']]/DHS_sample[['male_over_prop']]) /2
  ## Sample incidence from Lonnroth data, assuming independence of incidence estimates
  bmi_inc_est_sample <- rep(0,27)  
  for (x in 1:27){
    bmi_inc_est_lnorm      <- suppressWarnings(suppressMessages(get.lnorm.par(p=c(0.025,0.5,0.975),q=c(bmi_inc_rel[x,"inc_lo"],bmi_inc_rel[x,"inc"],bmi_inc_rel[x,"inc_hi"]),show.output=FALSE,plot=FALSE)))
    bmi_inc_est_sample[x]  <- rlnorm(1, meanlog = bmi_inc_est_lnorm[1], sdlog = bmi_inc_est_lnorm[2])
  }
  bmi_inc_sample <- data.frame(author=bmi_inc_rel$author, bmi=bmi_inc_rel$bmi, inc=bmi_inc_est_sample)

  ## 4. Fitting linear model =====================================================
  ## Dependent variable = incidence, independent variable = BMI, categorical variable = study
  ## Incidence = 10^(bmi_intercept + bmi_slope*bmi + study_coeff)
  model_bmi     <- lm(log10(inc) ~ bmi + author, data = bmi_inc_sample)
  bmi_slope     <- coef(summary(model_bmi))["bmi","Estimate"]
  
  ## 5. Plot the results to sense-check ========================================
  ## For use outside of the sampling loop to visually check results
  # bmi_intercept <- coef(summary(model_bmi))["(Intercept)","Estimate"]
  # q <- ggplot(bmi_inc_rel, aes(x=bmi, y=inc, colour=author))+
  #  geom_point()+
  ## Include error bars for incidence if available (not available in the sampling)
  #  geom_errorbar(aes(ymin=inc_lo,ymax=inc_hi))+
  ## Add linear model results in for each data set, where Cegielski is the reference
  #  geom_abline(intercept = bmi_intercept,                                                        slope=bmi_slope, colour="#00AD52", size=1)+
  #  geom_abline(intercept = bmi_intercept + coef(summary(model_bmi))["authoredwards","Estimate"], slope=bmi_slope, colour="#F62F31", size=1)+
  #  geom_abline(intercept = bmi_intercept + coef(summary(model_bmi))["authorhemila","Estimate"],  slope=bmi_slope, colour="#943094", size=1)+
  #  geom_abline(intercept = bmi_intercept + coef(summary(model_bmi))["authorleung","Estimate"],   slope=bmi_slope, colour="#595A5B", size=1)+
  #  geom_abline(intercept = bmi_intercept + coef(summary(model_bmi))["authorpalmer","Estimate"],  slope=bmi_slope, colour="#33399C", size=1)+
  #  geom_abline(intercept = bmi_intercept + coef(summary(model_bmi))["authortverdal","Estimate"], slope=bmi_slope, colour="#942D20", size=1)+
  #  ## Axis titles and limits
  #  scale_x_continuous(name="BMI (kg/m^2)",limits=c(15,35),minor_breaks = seq(15, 35, 1))+
  #  scale_y_log10(name="TB incidence/100k/year",limits=c(1,1000),minor_breaks = seq(0, 1000, 10))+
  #  scale_colour_manual(values = c("#00AD52", "#F62F31", "#943094", "#595A5B", "#33399C", "#942D20"))
  # ggsave("BMIrel.png",device="png",width=30,height=30,units=c("cm"))
  
  ## 6. Estimating relative risk by BMI ==========================================
  ## We use known incidence, BMI distribution and our linear model to estimate incidence by BMI
  ## We know that                                   inc_tot    = inc_modthin*prop_modthin + inc_mildthin*prop_mildthin + inc_normal*prop_normal + inc_over*prop_over
  ## We can use our linear model for two different known y (incidence) and known x (BMI) at those points, solving for the (unknown) intercept with (x1,y1) and substituting back into the equation for (x2,y2)
  ## Therefore:                                     inc_over   = 10^(log(inc_normal) + slope*(bmi_over-bmi_normal))
  ## So:                                            inc_tot    = 10^(log(inc_normal)+slope*(bmi_modthin-bmi_normal))*prop_modthin + 10^(log(inc_normal)+slope*(bmi_mildthin-bmi_normal))*prop_mildthin + 10^(log(inc_normal)+slope*(bmi_over-bmi_normal))*prop_over + inc_normal*prop_normal
  ## All of these are known except inc_normal, so inc_normal   = inc_tot/(prop_modthin*10^(slope*(bmi_modthin-bmi_normal)) +prop_mildthin*10^(slope*(bmi_mildthin-bmi_normal)) + prop_over*10^(slope*bmi_over-bmi_normal) + prop_normal)
  ## So for a given year's sampled incidence and BMI,
  inc_normal     <- inc_est_sample / (DHS_sample$modthin_prop*10^(bmi_slope*(DHS_sample$modthin_avg-DHS_sample$normal_avg)) + 
                                     DHS_sample$mildthin_prop*10^(bmi_slope*(DHS_sample$mildthin_avg-DHS_sample$normal_avg)) + 
                                     DHS_sample$over_prop*10^(bmi_slope*(DHS_sample$over_avg-DHS_sample$normal_avg))  +
                                     DHS_sample$normal_prop)
  inc_modthin    <- 10^(log10(inc_normal)+bmi_slope*(DHS_sample$modthin_avg-DHS_sample$normal_avg))
  inc_mildthin   <- 10^(log10(inc_normal)+bmi_slope*(DHS_sample$mildthin_avg-DHS_sample$normal_avg))
  inc_over       <- 10^(log10(inc_normal)+bmi_slope*(DHS_sample$over_avg-DHS_sample$normal_avg))
  ## Calculating incidence rate ratios compared to the overweight/obese category
  rr_modthin[y]  <- inc_modthin  / inc_normal
  rr_mildthin[y] <- inc_mildthin / inc_normal
  rr_over[y]     <- inc_over     / inc_normal
}
## Distribution of results
boxplot(rr_modthin,rr_mildthin,rr_over)
hist(rr_modthin)
hist(rr_mildthin)
hist(rr_over)

ci_modthin  <- mean(rr_modthin)  + qt( c(0.05, 0.95), sample_sz - 1) * sd(rr_modthin) 
ci_mildthin <- mean(rr_mildthin) + qt( c(0.05, 0.95), sample_sz - 1) * sd(rr_mildthin)
ci_over     <- mean(rr_over)     + qt( c(0.05, 0.95), sample_sz - 1) * sd(rr_over)

ci_tb_modthin  <- mean(tb_modthin)  + qt( c(0.05, 0.95), sample_sz - 1) * sd(tb_modthin) 
ci_tb_mildthin <- mean(tb_mildthin) + qt( c(0.05, 0.95), sample_sz - 1) * sd(tb_mildthin)
ci_tb_normal   <- mean(tb_normal) + qt( c(0.05, 0.95), sample_sz - 1) * sd(tb_normal)
ci_tb_over     <- mean(tb_over)     + qt( c(0.05, 0.95), sample_sz - 1) * sd(tb_over)

## 7. Trend in proportion of mod vs mild =======================================
## We use BMI and F:M estimates from statcompiler for 2005, 2015, 2019 https://www.statcompiler.com/en/
## Use this to fit a linear model for changing proportions over time
ModMild <- data.frame(year = c(2005, 2015, 2019),
                      MF = c(1.072, 1.069, 1.068),
                      F_mild = c(19.7, 13.3, 11.0),
                      F_mod = c(15.8, 9.6, 7.7),
                      M_mild = c(20.1, 12.0, 9.3),
                      M_mod = c(13.6, 7.8, 6.3))
ModMild$ratio <- (ModMild$F_mod/(ModMild$F_mod+ModMild$F_mild)+ModMild$MF*ModMild$M_mod/(ModMild$M_mod+ModMild$M_mild))/2
lm(ModMild$ratio ~ ModMild$year)











  