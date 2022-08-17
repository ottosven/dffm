# dffm


<!-- badges: start -->


<!-- badges: end -->


The goal of dffm is to provide functionality to apply the methods
developed in the paper “[Approximate Factor Models for Functional Time Series](https://arxiv.org/abs/2201.02532)” by [Sven
Otto](https://www.svenotto.com) and [Nazarii Salish](https://sites.google.com/site/nazariisalish/home).


## Preprint: 


https://arxiv.org/abs/2201.02532


## Installation

You can install the package using the following command:

``` r
library(remotes)
install_github("ottosven/dffm")
```

## Example
``` r
library(dffm)
# ---- data ---- #
data = load.fed()
# ---- preliminary --- #
fpca.obj = fpca.preprocess(data = data, method = "splines")
dffm.3Dplot(fpca.obj, domainlab = "maturities", outputlab = "yields (in percent)", 
            main = "3D plot of FED yields from 2001 to 2021")
plot(fpca.obj$eigenvalues, type = "b", main = "Screeplot", ylab = "eigenvalues")
dffm.criterion(fpcaobj = fpca.obj)
# ---- fitting yieldcurve ---- #
dffm.obj = dffm(fpcaobj = fpca.obj, K = 4, p = 1)
plot(dffm.obj$fittedcurve.workgrid[246,], type = "l", xlab = "maturities", ylab = "yields (in percent)", 
     main = "Yieldcurve of december 2021")
# ---- predicting the yieldcurve ---- #
predicted.dffm = dffm.forecast(dffmobj = dffm.obj, h = 10)
plot(predicted.dffm$predcurves.workgrid[10,], type = "l", xlab = "maturities", ylab = "yields (in percent)", main = "Predicted yieldcurve of october 2022")
```
