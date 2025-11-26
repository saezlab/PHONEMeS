---
layout: page
title: Installation
---

# PHONEMeS
For the installation of the most recent implementation of PHONEMeS the user can do this directly from GitHub by using the [devtools](https://www.r-project.org/nosvn/pandoc/devtools.html) package:

```R
# Install PHONEMeS from Github using devtools
install.packages('devtools') # in case devtools hasn't been installed
devtools::install_github('saezlab/PHONEMeS')
```

# PHONEMeS-stoch (not recommended)

For the installation of the old stochastic implementation (PHONEMeS-stoch), download the `tar`-file of the package (can be found [here](https://github.com/saezlab/PHONEMeS/tree/master/Package)). Open `R` and type:

```R
install.packages("PHONEMeS_0.2.7.tar.gz", repos=NULL)
```
