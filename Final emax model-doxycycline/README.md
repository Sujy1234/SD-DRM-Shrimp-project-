# Doxycycline pharmacodynamic parameter estimation

This folder contains the data and R code used to estimate the doxycycline pharmacodynamic parameters E<sub>max</sub> and EC<sub>50</sub> from published time-kill data.

The paired empirical bootstrap distributions of E<sub>max</sub> and EC<sub>50</sub> are used as inputs in the simple death dose-response model (SD-DRM).

## Analysis file

- `Emax code.R`  
  Main analysis script. It calculates concentration-specific net growth or killing rates from the time-kill data and fits an E<sub>max</sub> pharmacodynamic model with the Hill coefficient fixed at 1. A residual bootstrap with 1,000 iterations is used to quantify uncertainty in E<sub>max</sub> and EC<sub>50</sub>.

## Input data

- `Antibiotic data with control.csv`  
  Extracted time-kill data for doxycycline and the untreated control. The analysis uses observations collected up to 6 hours.

## Output files

- `Doxycycline_SD_DRM_parameter_summary.csv`  
  Point estimates and 95% bootstrap confidence intervals for E<sub>max</sub> and EC<sub>50</sub>.

- `Doxycycline_SD_DRM_empirical_bootstrap.csv`  
  Paired bootstrap estimates of E<sub>max</sub> and EC<sub>50</sub>. The paired rows are retained for use in the SD-DRM so that the empirical relationship between the two parameters is preserved.

## Supporting documents

- `Paper with time kill data.pdf`  
  Published study from which the doxycycline time-kill data were obtained.

- `Fitting time kill data into PD model.pdf`  
  Supporting document describing the pharmacodynamic modelling approach.

## Running the analysis

Place all files in the same working directory and run:

```r
source("Emax code.R")
