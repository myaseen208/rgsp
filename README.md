# rgsp: Repetitive Group Sampling Plan Based on Cpk

## Introduction
The application of statistical quality control tools in industry are gaining an increasing importance in the quality conscious modern society. Many new methods are developed in recent years but not commonly used in industry because of their non-availability in different commercial statistical softwares. An R package for repetitive group sampling plan proposed by Aslam et al. (2012) (http://dx.doi.org/10.1080/00949655.2012.663374) was developed and named as  ``rgsp``. The proposed sampling plan is based on process capability index Cpk, when the quality characteristic follows a normal distribution with unknown mean and variance. This paper presents methodology to calculate the minimum average sample number and average sample size in repetitive sampling for both symmetric and asymmetric cases by using ``rgsp`` package. The functions ``rgsp_asym1``, ``rgsp_asym2`` and ``rgsp_asym`` are used to estimate P1, P2, n, ka, kr  and ASN for both symmetric and asymmetric cases. 


## Authors
1. Muhammad Yaseen (myaseen208@gmail.com)
2. Muhammad Aslam (aslam_ravian@hotmail.com)
3. Sami Ullah (samiullahuos@gmail.com)
4. Muhammad Kashif (mkashif@uaf.edu.pk)


## Installation
Use **remotes** to install the development version from Github:

```{r}
if(!require("remotes")) install.packages("remotes")
remotes::install_github("myaseen208/rgsp")
```

## License
This package is free and open source software, licensed under GPL.
