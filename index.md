---
layout: default
title: Home
---


# PHONEMeS

## Overview

PHONEMeS (**PHO**sphorylation **NE**tworks for **M**ass **S**pectrometry) is a method to model signaling networks based on untargeted phosphoproteomics mass spectrometry data and kinase/phosphatase-substrate interactions. This is a tool currently being developed by the [saezlab](https://saezlab.org/) members. The tool has been implemented as an R-package which contains the accompanying R functions as well as several examples about how to run a PHONEMeS analysis. The method was originally introduced from [Terfve et.al. 2015](https://www.nature.com/articles/ncomms9033) while major improvements and updates have been provided in [Gjerga et.al. 2021](https://pubs.acs.org/doi/abs/10.1021/acs.jproteome.0c00958).

According to the pipeline shown in (*Figure 1*), cells are treated by a combination of kinase inhibitors which target specific intracellular proteins or receptor activators. After treatment, samples of tissue undergo a discovery mass-spectrometry based analysis where phospho-proteomic changes are quantified for each peptide at each condition. Stringent
statistical analysis of the mass-spectrometry data is performed in order to identify the effects that perturbation treatments have over each phospho-site. These analyses typically consist of a differential expression analysis between treated and control samples from where regulation levels (i.e. log fold changes). Robustness of such regulation can be also estimated through statistical testing. A score is assigned to each of the measured phospho-proteins at each experimental condition based on the regulation levels and the estimated statistics. Results from differential expression analysis are then combined with prior knowledge about signalling. The prior knowledge typically consists of kinase-substrate
interactions but it can as well be complemented with directed protein-interactions. Training of the prior knowledge network to data has been recently formulated as an integer linear programming program ([Gjerga et.al. 2021](https://pubs.acs.org/doi/abs/10.1021/acs.jproteome.0c00958)) enhancing model optimisation speed by orders of magnitude compared to the previous genetic algorithm implementation. After training, the most likely linear paths are obtained and which are linking our perturbation targets to the regulated phosphosites identified by the analysis of the phosphoproteomics datasets (example, *Figure 2*). Our contextualised network solutions can then be used to generate testable hypotheses and propose possible therapeutic treatments.

<img src="/PHONEMeS/public/phonemes_pipeline.png" alt="Example network">
<center><i>Figure1: PHONEMeS pipeline.</i></center> 

Information about the various PHONEMeS variants implemented and usage details about PHONEMeS-ILP can be found [here](https://saezlab.github.io/PHONEMeS/2_usage/).

<img src="/PHONEMeS/public/schaefer_network.jpg" alt="Example network">

<center><i>Figure 2: Time-resolved network model from PHONEMeS for UACC257 cells, figure from [Schafer et.al. 2019](https://www.embopress.org/doi/full/10.15252/msb.20198828) (licensed under the Creative Commons Attribution CC BY 4.0 ). EDNRB (green diamond) was connected to its target phosphorylation sites (red hexagons) through intermediary kinases or G proteins (blue circles) in a time‐resolved variant of the PHONEMeS approach. Central kinases, which were also identified by kinase activation prediction, are shown as intermediary kinases with dark blue shading. Edge thickness corresponds to weights, which were assigned by downsampling the network 100 times.</i></center>

## References

If you use any of the ILP implementations of PHONEMeS (PHONEMeS-ILP/PHONEMeS-mult/PHONEMeS-UD), please cite:

> Gjerga E., Dugourd A., Tobalina L., Sousa A., Saez-Rodriguez J. (2021). [PHONEMeS: Efficient Modeling of Signaling Networks Derived from Large-Scale Mass-Spectrometry Data](https://pubs.acs.org/doi/full/10.1021/acs.jproteome.0c00958) _J Prot Res_, 2021, DOI: 10.1021/acs.jproteome.0c00958.

If you use the original implementation of PHONEMeS or if you would like to make mention of the tool, please additionally cite:

> Terfve C. D. A., Wilkes E. H., Casado P., Cutillas P. R., and Saez-Rodriguez J. (2015). [Large-scale models of signal propagation in human cells derived from discovery phosphoproteomic data.](http://www.nature.com/articles/ncomms9033) _Nature Communications_, **6**:8033.


### Other References
 
+ Application of PHONEMeS-ILP over time-resolved data:

  > Schäfer A., Gjerga E., Welford R.W.D., Renz I., Lehembre F., Groenen P.M.A., Saaez-Rodriguez J., Aebersold R., and Gstaiger M. (August 2019). [Elucidating essential kinases of endothelin signalling by logic modelling of phosphoproteomics data](https://www.embopress.org/doi/abs/10.15252/msb.20198828) _Molecular Systems Biology_, **Volume 15, Issue 8**. Codes for this project can be found [here](https://github.com/saezlab/EDN_phospho).

+ Application of PHONEMeS-ILP to analyze the effect of hypertensive insults on kidneys:

  > Rinschen M.M., Palygin O., Guijas C., Palermo A., Palacio-Escat N., Domingo-Almenara X., Montenegro-Burke R., Saez-Rodriguez J., Staruschenko A. and Siuzdak G. (December 2019). [Metabolic rewiring of the hypertensive kidney](https://stke.sciencemag.org/content/12/611/eaax9760) _Science Signalling_, **Vol. 12, Issue 611, eaax9760**
 
 + Application of PHONEMeS-ILP over colon cancer data:
  > Gjerga E., Dugourd A., Tobalina L., Sousa A., Saez-Rodriguez J. (2021). [PHONEMeS: Efficient Modeling of Signaling Networks Derived from Large-Scale Mass-Spectrometry Data](https://pubs.acs.org/doi/full/10.1021/acs.jproteome.0c00958) _J Prot Res_, 2021, DOI: 10.1021/acs.jproteome.0c00958.
