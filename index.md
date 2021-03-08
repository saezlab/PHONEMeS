---
layout: default
title: Home
---


# PHONEMeS

## Overview

PHONEMeS (**PHO**sphorylation **NE**tworks for **M**ass **S**pectrometry) is a method to model signaling networks based on untargeted phosphoproteomics mass spectrometry data and kinase/phosphatase-substrate interactions.

This package contains the R package and accompanying scripts that implement the method as well as several examples of how to run a PHONEMeS analysis.

The input for PHONEMeS consists of phosphoproteomic data after treatment with kinase inhibitors. Gaussian mixture modeling is then used to find phosphosites that exhibit a naturally Boolean behaviour with two populations, representing a control and a perturbed state. The data are mapped unto a kinase/phosphatase-substrate network taken from several dedicated databases. Optionally the users can also rely on standard differential analysis results of phospho-proteomics data as an input. PHONEMeS then optimizes the network and extracts possible paths connecting inhibited kinases and perturbed phosphosites. In the [original implementation](https://www.nature.com/articles/ncomms9033) is performed by iteratively sampling edges from the background network, simulating the logic model, and finally evaluating the network by comparison to the data. In each iteration, sampling weights are corrected based on the edge frequencies in the best models. 

**Recently PHONEMeS has been has been reformulated as an Integer Linear Programming (PHONEMeS-ILP) which is orders of magnitude faster than the original implementation, it requires no cluster infrastructure to run the analysis and expands the scenarios that can be analyzed**

Usage details about PHONEMeS-ILP can be found [here](https://saezlab.github.io/PHONEMeS/2_usage/).

<img src="/PHONEMeS/public/network.png" alt="Example network">

<center><i>Example output of PHONEMeS analysis after MTOR inhibition</i></center>

## References

If you use PHONEMeS-ILP, please cite:

> Gjerga E., Dugourd A., Tobalina L., Sousa A., Saez-Rodriguez J. (2021). [PHONEMeS: Efficient Modeling of Signaling Networks Derived from Large-Scale Mass-Spectrometry Data](https://saezlab.org/) _J Prot Res_, 

If you use the original implementation of PHONEMeS or if you would like to make mention of the tool, please additionally cite:

> Terfve C. D. A., Wilkes E. H., Casado P., Cutillas P. R., and Saez-Rodriguez J. (2015). [Large-scale models of signal propagation in human cells derived from discovery phosphoproteomic data.](http://www.nature.com/articles/ncomms9033) _Nature Communications_, **6**:8033.


### Other References

+ Description of parts of the data:

  > Wilkes, E. H., Terfve, C., Gribben, J. G., Saez-Rodriguez, J., and Cutillas, P. R. (2015). [Empirical inference of circuitry and plasticity in a kinase signaling network.](http://www.pnas.org/content/112/25/7719.abstract) _Proceedings of the National Academy of Sciences of the United States of America_, **112**(25):7719–24.
 
+ Application of PHONEMeS-ILP over time-resolved data:

  > Schäfer A., Gjerga E., Welford R.W.D., Renz I., Lehembre F., Groenen P.M.A., Saaez-Rodriguez J., Aebersold R., and Gstaiger M. (August 2019). [Elucidating essential kinases of endothelin signalling by logic modelling of phosphoproteomics data](https://www.embopress.org/doi/abs/10.15252/msb.20198828) _Molecular Systems Biology_, **Volume 15, Issue 8**. Codes for this project can be found [here](https://github.com/saezlab/EDN_phospho).

+ Application of PHONEMeS-ILP to analyze the effect of hypertensive insults on kidneys:

  > Rinschen M.M., Palygin O., Guijas C., Palermo A., Palacio-Escat N., Domingo-Almenara X., Montenegro-Burke R., Saez-Rodriguez J., Staruschenko A. and Siuzdak G. (December 2019). [Metabolic rewiring of the hypertensive kidney](https://stke.sciencemag.org/content/12/611/eaax9760) _Science Signalling_, **Vol. 12, Issue 611, eaax9760**
 
 + Application of PHONEMeS-ILP over colon cancer data:
 > Gjerga E., Dugourd A., Tobalina L., Sousa A., Saez-Rodriguez J. (2021). [PHONEMeS: Efficient Modeling of Signaling Networks Derived from Large-Scale Mass-Spectrometry Data](https://saezlab.org/) _J Prot Res_, 
