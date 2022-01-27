# PHONEMeS

**PHONEMeS** (**PHO**sphorylation **NE**tworks for **M**ass **S**pectrometry) is a method to model signalling networks based on untargeted phosphoproteomics mass spectrometry data and kinase/phosphatase-substrate interactions. 
The newest version of PHONEMeS builds on `CARNIVAL`, and more specifically the `runCARNIVAL`wrapper function. As such, it requires the CARNIVAL package (Version 1.3) to be installed. Additionally, it uses `OmnipathR` to construct prior knowledge networks based on solely kinase-substrate interactions or the combination of kinase-substrate and protein-protein interactions.
The pipeline requires a vector of input nodes (perturbed kinases) and measured nodes (deregulated phosphosites) that are connected by an interactive version of 
IBM Cplex or CBC-COIN solver. The IBM ILOG Cplex is freely available through Academic Initiative [here](https://www.ibm.com/products/ilog-cplex-optimization-studio). The [CBC](https://projects.coin-or.org/Cbc) solver is open source and freely available
for any user. Alternatively for smaller cases, users can rely on the freely available 
[lpSolve R-package](https://cran.r-project.org/web/packages/lpSolve/index.html). 


### License

Distributed under the GNU GPLv2 License. See accompanying file [LICENSE.txt](https://github.com/saezlab/PHONEMeS/blob/master/LICENSE.txt) or copy at https://www.gnu.org/licenses/gpl-2.0.html.

### Installation

To install `PHONEMeS` please run:
```
devtools::install_github('saezlab/PHONEMeS')
```

### Usage

A short tutorial on how to run a PHONEMeS analysis using PHONEMeS v2.0.0 is provided
in the [vignette tutorial](https://github.com/saezlab/PHONEMeS/blob/master/vignettes/tutorial.Rmd).


### Prior versions

The code for the original PHONEMeS package (PHONEMeS v1.0.0) as described in Terfve et al. 2015 can be found in the releases.
For a guide how to run a PHONEMeS analysis using PHONEMeS v1.0.0, please refer to the [documentation](https://saezlab.github.io/PHONEMeS).


### References

[Terfve et al.](http://www.nature.com/articles/ncomms9033):

> Terfve, C. D. A., Wilkes, E. H., Casado, P., Cutillas, P. R., and Saez-Rodriguez, J. (2015). Large-scale models of signal propagation in human cells derived from discovery phosphoproteomic data. *Nature Communications*, 6:8033.

[Wilkes et al.](http://www.pnas.org/content/112/25/7719.abstract) (description of parts of the data)

> Wilkes, E. H., Terfve, C., Gribben, J. G., Saez-Rodriguez, J., and Cutillas, P. R. (2015). Empirical inference of circuitry and plasticity in a kinase signaling network. *Proceedings of the National Academy of Sciences of the United States of America,* 112(25):7719–24.

