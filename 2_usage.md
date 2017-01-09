---
layout: page
title: Usage
---

This documentation is based on the scripts and data files provided in the [PHONEMeS Main Example folder](https://github.com/saezlab/PHONEMeS/tree/master/Example/Example_MainData).

The documentation is based on PHONEMeS 0.2.3, run locally using R 3.1.0 (in RStudio) and on a cluster using R 2.15.2 and LSF 8. 

These scripts require the *[igraph](http://igraph.org/r/)* package (version 0.7.1) and the *[BioNet](https://www.bioconductor.org/packages/release/bioc/html/BioNet.html)* package (version 1.24.0). 

The result networks are visualized using *[Cytoscape](http://www.cytoscape.org/)* 2.8.0 (although we are working on a plugin for Cytoscape 3 to automatically import the contents of a typical PHONEMeS results folder and to produce annotated networks).

# PHONEMeS workflow

# I. Data Preparation

The untargeted phosphoproteomics data used for the identification of the signaling models in PHONEMeS is a data matrix containing the peak heights for each of the measured peptides accross each condition (treatments and control) and replicates and which maps to specific sites over a certain number of proteins. If the data is not normalized, we quantile normalize the *log10* of the raw intensity values and then use a linear model to estimate the effects of each drug over a control state by computing the log fold change between each point of these two conditions and estimating their significance using a *t-statistic* as computed by the function *eBayes* implemented by the *Bioconductor* package *[limma](http://www.bioconductor.org/packages/2.12/bioc/html/limma.html)*

As a second step, we fit a Gaussian Mixture Model (GMM) on each peptide accross each condition and select only those peptides whose distributions are better described with 2 components, hence showing a Boolean behaviour under each different condition (perturbed/non-perturbed). For that we can use the *Bioconductor* package *[mclust](http://www.stat.washington.edu/mclust/)* which runs an expectation-maximization method to fit the measurements into mixed Gaussian distributions. Furthermore, we exclude those cases in which the two density curves overlap by more than *10%* over a range covering the whole data in order to avoid cases where two conditions are found because of a tightly clustered set of measurements. Then for each point, we associate the score ![first equantion](http://latex.codecogs.com/gif.latex?S_%7Bi%2C%20j%7D) as the log ratio of the probabilty of belonging to either a control or perturbed state: ![second equation](http://latex.codecogs.com/gif.latex?S_%7Bi%2C%20j%7D%20%3D%20%5Clog10%28%20%5Cfrac%7B%20P_%7Bi%2C%20j%7D%20%28%20C_%7Bi%7D%20%29%7D%7BP_%7Bi%2C%20j%7D%20%28%20P_%7Bi%7D%20%29%7D%20%29)

Peptide *i* is called perturbed under condition *j* if ![third equation](http://latex.codecogs.com/gif.latex?S_%7Bi%2C%20j%7D%20%3C%20-0.5) and it is considered to be in the control state if ![fourth equation](http://latex.codecogs.com/gif.latex?S_%7Bi%2C%20j%7D%20%3E%200.5). Values in between *-0.5* and *0.5* are considered as undetermined.

**TODO_2**: Update this documentation with the new scripts from cluster

# II. Running PHONEMeS on a Cluster

### 1. Prepare the data for the cluster

Run `prepOptim.R` and produce `data4cluster_n.RData`. Here, `n` denotes the index of the independent optimization run.
`prepOptim.R` takes a network object and a data object (resulting from data normalization, summarization, and Gaussian mixture modeling) in order to produce the objects needed to run PHONEMeS on a cluster.

### 2. Move to cluster

Copy the following files to a cluster:

 * data4cluster_n.RData
 * processGx_n.R
 * runScriptGx_n.sh
 * scriptGxopt_50models_n.R
 * import_n.R

### 3. Is the package installed?

Make sure that PHONEMeS is installed on your cluster. If this is not the case, run the following in R:

```R
install.packages("PHONEMeS_0.2.3.tar.gz", repos=NULL)
```

### 4. Create results directory

Additionally, each independent optimization run requires a result folder with the index `n` of the optimization run.

```bash
mkdir Results_n
```

### 5. Run one independent optimization

This step should be repeated multiple times with different indices. The script `runScriptGx_n.sh` will submit the necessary jobs to the queue of your LSF system. Make sure to modify the script as necessary (e.g. correcting the queue name).

```bash
./runScriptGx_n.sh
```

### 6. Move results to your local machine

Copy back the following files back to your local machine in order to process, combine, and visualize the results of the PHONEMeS analysis:

 * pn_imported.RData
 * optim_n.pdf

# III. Post-Optimization Analysis

### Process each independent optimization

The R script `postOptim.R` will process the results of a single optimization, so that it should be run for each independent optimization run. The file may require minor modifications regarding file location.

### Combine multiple independent optimizations

Once `postOptim.R` has been run on all independent optimization and produced the different files `objects_pn.RData`, run the R script `comb_optim.R`. It will output the combined plots as well as the final resulting network (maximal input averaged frequencies across independent optimizations, maximal scoring paths by averaged frequencies, etc.).

# IV. Visualize the Result Networks

The resulting `.txt`-files (either individual ones from [single](#process-each-independent-optimization) or combined [optimizations](#combine-multiple-independent-optimizations)) can be imported as tables into cytoscape.

1. Start Cytoscape and select "From Network File...".

2. Select `nTagMaxIn_comb.txt` and set the columns `K.ID` as source and `S.cc` as target.

3. Import edge frequencies: 
  * "Import Table from File": Select `combOptim_EA.txt`. 
  * Select "Import Data as: Edge Table Column" and "Key Column for Network: SID". 
  * Select the columns `SID` as Key and the column `f50` as Attribute.

4. Import Nodes Attributes:
  * "Import Table from File": Select `AllNodes_nodesP_NA_pn.txt`. 
  * Select "Import Data as: Node Table Column".

5. Import visual properties:
  * "Import Styles": Select `PHONEMeS_vizmap.props` in the PHONEMeS repository ([here](https://github.com/saezlab/PHONEMeS/tree/master/Example)).
  * Finally, select the PHONEMeS style in the "Style" tab.
