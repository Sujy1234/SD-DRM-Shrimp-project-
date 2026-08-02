# SD-DRM shrimp project

This repository contains the code and data used to develop and apply a simple death dose-response model (SD-DRM) for antibiotic-resistant *Vibrio parahaemolyticus* associated with imported shrimp.

## Repository structure

- `Final emax model-doxycycline`  
  Contains the doxycycline time-kill data and R code used to estimate the pharmacodynamic parameters \(E_{\max}\) and \(EC_{50}\). Residual bootstrap analysis is used to generate paired empirical estimates of these parameters for use in the SD-DRM.

- `SD-DRM(SHRIMP)`  
  Contains the main shrimp SD-DRM analysis and sensitivity analysis. The model estimates infection risks for raw shrimp consumption, undercooked shrimp consumption, and hand-mediated cross-contamination, and evaluates how the resistant fraction and residual doxycycline concentration influence infection risk and predicted treatability.

The doxycycline pharmacodynamic parameters estimated in the first folder are used as inputs in the shrimp SD-DRM analysis.
