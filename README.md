# Reproducibility for Sung and Hung (2023) Efficient Calibration for Imperfect Epidemic Models with Applications to the Analysis of COVID-19
This folder consists of the data and R code for the manuscript [Efficient Calibration for Imperfect Epidemic Models with Applications to the Analysis of COVID-19](https://arxiv.org/abs/2009.12523) by Sung and Hung (2023). 

* The code folder reproduces the results in the Sections 4 and 5 of the manuscript. 
  * The required R packages include `hetGP`, `numDeriv`, `lhs`, `SimInf`, `deSolve`, `covid19.analytics`, `ggplot2`, `snowfall`, `doParallel`, `foreach`, `plyr`, `gelnet`, `MRFA`, `grplasso`, `scales`, and `GPfit`
  * `numericalstudy_1d.R` reproduces Figure 1.
  * `numericalstudy_3d.R` reproduces Figure 2.
  * `COVID19_analysis.R` reproduces Figures 3-7. WARNING: this code could take more than 24 hours to run, depending on the computer resourses.

