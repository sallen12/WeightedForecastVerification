# WeightedForecastVerification

This repository contains R code to reproduce the results presented in the paper  

> Allen, S., Bhend, J., Martius, O. and Ziegel, J. (2023). 
> Weighted verification tools to evaluate univariate and multivariate probabilistic forecasts for high-impact weather events
> Weather and Forecasting. 38, 499–516.
> [https://doi.org/10.1175/WAF-D-22-0161.1](https://doi.org/10.1175/WAF-D-22-0161.1)

## Weighted verification tools

This paper discusses how user-specified weight functions can be incorporated into forecast evaluation, allowing forecasters to assess their predictions whilst emphasising particular outcomes.

Probabilistic forecasts are typically assessed with respect to their accuracy and their calibration. 

### Weighted scoring rules

Scoring rules provide a measure of forecast accuracy, and allow competing forecasters to be ranked and compared objectively.
Weighted scoring rules have become a well-established tool with which to emphasise particular outcomes within conventional scoring rules.
Popular weighted scoring rules are available in the [scoringRules package](https://github.com/FK83/scoringRules), but only for forecasts in 
the form of predictive samples. In this repository, a collection of weighted scoring rules is provided for familiar univariate parametric distributions.
Currently, this is available for the normal, logistic, and Student's t distributions, though more could be added in the future.

### Weighted PIT histograms

Forecast calibration is typically assessed using rank and probability integral transform (PIT) histograms. 
The above paper illustrates how the theory underlying weighted scoring rules can be extended to checks for forecast calibration.
This repository contains code to calculate conditional PIT values, and to plot rank histograms, PIT histograms, and PIT reliability diagrams.

The vignette thoroughly documents the usage of these weighted verification tools, and reproduces the results in the above paper.

## Data

The data used in this study is not publicly available, so the repository instead contains a 'noisy' data set that seeks to artificially mimic the features of the actual data. Hence, the results presented in the vignette do not correspond exactly to those presented in the paper. Further details regarding the data can be found in the paper and the vignette.

## Installation and development

This package has not been submitted to CRAN, and can therefore be installed in R using devtools
```r
# install.packages("devtools")
library(devtools)
install_github("sallen12/WeightedForecastVerification")
```
The package is still in active development, and the vignette lists several possible extensions that could be implemented. Comments, suggestions, and input are more than welcome.
