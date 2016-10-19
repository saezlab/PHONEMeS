# PHONEMeS

**PHONEMeS** (**PHO**sphorylation **NE**tworks for **M**ass **S**pectrometry) is a method to model signalling networks based on untargeted phosphoproteomics mass spectrometry data and kinase/phosphatase-substrate interactions. 
Please see [Terfve et al.](http://www.nature.com/articles/ncomms9033) for an explanation of the methodolgy and as an example for how to run a PHONEMeS analysis.

This repository contains the [R package](https://github.com/saezlab/PHONEMeS/tree/master/Package) and accompanying scripts that implement the method.

### License

Distributed under the GNU GPLv2 License. See accompanying file [LICENSE.txt](https://github.com/saezlab/PHONEMeS/blob/master/LICENSE.txt) or copy at https://www.gnu.org/licenses/gpl-2.0.html.

### Installation

For installation, download the tar file of the package and type in R:

```R
install.packages("PHONEMeS_0.2.3.tar.gz", repos=NULL)
```

### Usage

For a guide how to run a PHONEMeS analysis, please refer to the [documentation](https://github.com/saezlab/PHONEMeS/blob/master/Documentation/how_to.md).

### References

[Terfve et al.](http://www.nature.com/articles/ncomms9033):

> Terfve, C. D. A., Wilkes, E. H., Casado, P., Cutillas, P. R., and Saez-Rodriguez, J. (2015). Large-scale models of signal propagation in human cells derived from discovery phosphoproteomic data. *Nature Communications*, 6:8033.

[Wilkes et al.](http://www.pnas.org/content/112/25/7719.abstract) (description of parts of the data)

> Wilkes, E. H., Terfve, C., Gribben, J. G., Saez-Rodriguez, J., and Cutillas, P. R. (2015). Empirical inference of circuitry and plasticity in a kinase signaling network. *Proceedings of the National Academy of Sciences of the United States of America,* 112(25):7719–24.

