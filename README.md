# DataRemix
DataRemix: a universal data transformation for optimal inference from gene expression datasets


R library dependencies

```
MASS
ROCR
caTools
knitr(only required by DataRemix_display())
```

Extra library dependencies
```
Intel MKL or other BLAS/LAPACK library is highly recommended. For some large matrices, 
it can achieve about 70-fold acceleration.
```

To install the R package
```
library(devtools)
install_github("wgmao/DataRemix")
```

Check out the vignette [here](vignettes/vignette.pdf)

