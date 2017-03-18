---
layout: page
title: Examples
---

In this section we describe an example on how to run PHONEMeS on the original [data](https://github.com/saezlab/PHONEMeS/tree/master/Example/Example_MainData/data) from the original [publication](http://www.nature.com/articles/ncomms9033).


First, we run the ```prePHONEMeS_rawData2PHONEMeS.R``` script from which we naormalize our raw data ([MSdata.RData](https://github.com/saezlab/PHONEMeS/blob/master/Example/Example_MainData/data/MSdata.RData)), perform some statistical analysis of the normalized data by applying a linear model (see [Limma](http://bioconductor.org/packages/release/bioc/html/limma.html) package) thus obtaining the significance of the fold changes between the perturbed and control measurements. Then a Gaussian Mixture Model is applied (see [mclust](https://cran.r-project.org/web/packages/mclust/index.html) R package) in order to identify those peptides which are best fitted with a certain number of components/clusters. In addition we remove those peptides where the control is NA since we cannot assign a control cluster on this case. Only those measurements which are best described with two comoponents are taken into consideration. A probability is given to each of the peptides on each of the condition based on the resulting GMM parameters which tells how likely is this measurememnt to belong to the control or the perturbed state. The output of this script will be a data matrix ```GMM.res``` where we identify categories of data points after linear model estimation of fold changes vs control and Gaussian mixture modelling by associateing to the fold changes for drug vs control an adjusted p-value which is estimated with the linear model.

Run ```preOptim.R``` which takes as an input the data matrix output from ```prePHONEMeS_rawData2PHONEMeS.R```and a frame of data containing all possible K/P-Substrate interactions that we can retrieve from large data-bases such as [Omnipath](http://omnipathdb.org) (see original publication of Omnipath [here](http://www.nature.com/nmeth/journal/v13/n12/full/nmeth.4077.html)). After defining our target kinases(MTOR on this case), the script will choose the drug treatments matching to the drug targets and match to what is present in the background network. This will allow the building of a background prior knowledge network containing all possible interactions connecting our drug targets with the measurements we get. The ```preOptim.R``` will then create the objects needed to run PHONEMeS on a cluster and the optimization parameters such as:
* the size penalty factor *sizeP*
* number of models sampled in each script *nG1*
* a logical indicating whether should the direct interactions between drug targets and data sites be added to every model *cstart*
* another logical indicating whether we want integrators to be sampled the same way intermediates are *intgAsintm*
* number of parallel scripts that are run at each generation *nScripts*
* a logical indicating whether should absolute tolerance be used *absTol*
* the tolerance *tol*
* the maximum number of times that an edge is copied through at each generation *cap*
* index of the optimisation (for naming files and folders) *resN*
The number of the independent optimization runs on this example as we can see from the script is 6.
The output of the ```preOptim.R``` will be ```data4cluster_n.RData``` where ```n``` denotes the index of the independent optimization run.

Before running the analysis on the cluster, we have to make sure that the [PHONEMeS package](https://github.com/saezlab/PHONEMeS/blob/master/Package/PHONEMeS_0.2.7.tar.gz) is already installed on the cluster. After we make sure that PHONEMeS has been successfully installed, we import on the cluster the outputs from the previous step ```data4cluster_n.RData``` and the other following scripts (**important note:** in the cluster we create different directories for each independent run through ```mkdir``` command):
* ```processGx.R``` which combines all the results (model and scores) for each independent optimization runs.
* ```runScriptGx.sh``` for running the parallel processes
* ```scriptGxopt_50models.R``` which optimizes all the models sampled.
* ```import.R``` for locally importing the outputs from the analysis on the cluster.
After running the bash script```./runScriptGx.sh n```, on each of the folders we created we will have our outputs which we need to import locally.

After every result has been imported successfully, we process our results with ```postOptim.R``` which will assign attributes to each of the nodes and edges of each of the independent runs and save everything on one single folder (```objects_p3.RData```). Through ```comb_optim.R``` we will then combine all the ```objects_pn.RData``` in a sngle final resulting network.

As for the visualization of the resulting network, we can use Cytoscape. For more details check the description in the [Usage](https://saezlab.github.io/PHONEMeS/2_usage/) section.


