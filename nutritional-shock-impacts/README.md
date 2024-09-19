## Nutritional_Shock_impacts

This is a description for generating the output for *"Modelling the effect of nutritional shocks from the war in Ukraine on tuberculosis in India".*

***System and software requirements***

This project uses R version 4.3.2.

This code has been tested on macOS: Sonoma (14.4.1).

***R packages required***

```{r}
arrow
assertthat
cowplot
data.table
deSolve
digest
fst
getopt
ggplot2
here
hmer
lhs
log4r
log4r
logger
lubridate
Matrix
minpack.lm
patchwork
renv
stringi
tools
xml2
```

***Demo and instructions for use***

Country-specific files to run and calibrate the model are located in `TBmodel/countries/INDu/parameters` and country-specific demographic, mortality, and BMI distribution data are located in `TBmodel/countries/INDu/data`.

1.  Model calibration

    `TBmodel/0_AutoEmulate.R`

    -   Script that runs model calibration using history matching with emulation (`hmer`) to generate parameter sets that fit all targets

    -   Estimated length of time: 1-2 weeks to find 1000 fully fitting parameter sets. Given the length of time required to find parameter sets, this step was completed on a high-performance computing cluster.

2.  Model simulation

    `TBmodel/1_GenEpiOutput.R`

    -   Using the fully fitting parameter sets generated from calibration, this script will run the shock scenarios and generate the output

    -   Estimated length of time: 1 day to generate raw output

3.  Post-process model output

    `TBmodel/2_GenEpiUncertainty.R`

    -   Using the raw output generated from `2_GenEpiUncertainty.R`, this script will produce uncertainty intervals

    -   Estimated length of time: \<1 minute

4.  Plots

    `TBmodel/3_GenEpiPlots.R`

    -   Using the post processed output with uncertainty intervals from the previous step, the script will generate plots

    -   Estimated length of time: \<1 minute
