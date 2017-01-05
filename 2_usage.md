---
layout: page
title: Usage
---

This documentation is based on the scripts and data files provided in the [PHONEMeS Main Example folder](https://github.com/saezlab/PHONEMeS/tree/master/Example/Example_MainData).

The documentation is based on PHONEMeS 0.2.3, run locally using R 3.1.0 (in RStudio) and on a cluster using R 2.15.2 and LSF 8. 

These scripts require the *[igraph](http://igraph.org/r/)* package (version 0.7.1) and the *[BioNet](https://www.bioconductor.org/packages/release/bioc/html/BioNet.html)* package (version 1.24.0). 

The result networks are visualized using *[Cytoscape](http://www.cytoscape.org/)* 2.8.0 (although we are working on a plugin for Cytoscape 3 to automatically import the contents of a typical PHONEMeS results folder and to produce annotated networks).

# PHONEMeS workflow

## Data Preparation

## Running PHONEMeS on a Cluster

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

## Post-Optimization Analysis

### Process each independent optimization

The R script `postOptim.R` will process the results of a single optimization, so that it should be run for each independent optimization run. The file may require minor modifications regarding file location.

### Combine multiple independent optimizations

Once `postOptim.R` has been run on all independent optimization and produced the different files `objects_pn.RData`, run the R script `comb_optim.R`. It will output the combined plots as well as the final resulting network (maximal input averaged frequencies across independent optimizations, maximal scoring paths by averaged frequencies, etc.).

## Visualize the Result Networks

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
