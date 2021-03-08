---
layout: page
title: Installation
---

# PHONEMeS-ILP
For the installation of the most recent ILP implementation of PHONEMeS the user can do this directly from GitHub by using the [devtools](https://www.r-project.org/nosvn/pandoc/devtools.html) package:

```R
# Install PHONEMeS from Github using devtools
install.packages('devtools') # in case devtools hasn't been installed
library(devtools)
install_github('saezlab/PHONEMeS-ILP')
```

Otherwise users can install PHONEMeS directly from the source after downloading the source (tar file) and typing in ```R``` command line the following:

```R
# or download the source file from GitHub and install from source
install.packages('path_to_extracted_PHONEMeS_directory', repos = NULL, type="source")
```

For PHONEMeS-ILP depedencies, please read carefully the documentation of the [PHONEMeS-ILP package](https://github.com/saezlab/PHONEMeS-ILP).

# PHONEMeS-stoch

For the installation of the old stochastic implementation (PHONEMeS-stoch), download the `tar`-file of the package (can be found [here](https://github.com/saezlab/PHONEMeS/tree/master/Package)). Open `R` and type:

```R
install.packages("PHONEMeS_0.2.7.tar.gz", repos=NULL)
```
