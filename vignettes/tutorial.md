# Introduction

PHONEMeS is a method to model signaling networks based on untargeted
phosphoproteomics mass spectrometry data and
kinase/phosphatase-substrate interactions. It identifies a subset of
interactions from a prior knowledge network that represent potential
regulatory signaling pathways linking known or potentially deregulated
kinases to downstream phorsphorylation sites.

For the extraction of perturbed signaling networks using PHONEMeS, the
following inputs are required: 1. named vector of deregulated
phosphorylation sites 2. named vector of deregulated kinases 3. path to
a solver

The potentially deregulated phosphorylation sites and kinases will then
be connected based on a prior knowledge network. This prior knowledge
network can either be taken from within the PHONEMeS package or manually
added by the user. Within PHONEMeS there are two different prior
knowledge networks. One is based on the combination of kinase-substrate
interactions and protein-protein interactions (phonemesPKN). This
network is directed and signed where phosphorylation sites are always
end points and not connected to their downstream protein, as the effect
of a phosphorylation site is often unknown. This is the default network
that is used in run_phonemes. The other network is solely based on
kinase-substrate interactions and is only directed but not signed
(phonemesKSN). This network was also used in PHONEMeS 1.0 and is the
default network in run_vanilla_phonemes.

Here you can find a brief PHONEMeS tutorial. It provides an example for
an analysis with PHONEMeS starting from the results of an differential
analysis of phosphorylation sites, in this particular case a top table
results object from limma. This should allow users to prepare the input
objects for PHONEMeS starting from their own data.

# Getting set up

First, we load the PHONEMeS library. To use PHONEMeS, in particular the
run_phonemes function, an interactive version of IBM Cplex or CBC-COIN
solver is required as the network optimizer. The IBM ILOG Cplex is
freely available through Academic Initiative
[here](https://www.ibm.com/products/ilog-cplex-optimization-studio). The
[CBC](https://projects.coin-or.org/Cbc) solver is open source and freely
available for any user. Alternatively for smaller cases, users can rely
on the freely available [lpSolve
R-package](https://cran.r-project.org/web/packages/lpSolve/index.html).

Additionally, we load the
[decoupleR](https://github.com/saezlab/decoupleR) library, a recent
development from the Saez-Rodriguez group which offers 12 different
statistical approaches to perform functional analysis of omics data.

``` r
library(PHONEMeS)
library(decoupleR)
library(tidyverse)
#> ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──
#> ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
#> ✓ tibble  3.1.6     ✓ dplyr   1.0.8
#> ✓ tidyr   1.2.0     ✓ stringr 1.4.0
#> ✓ readr   2.1.2     ✓ forcats 0.5.1
#> Warning: package 'tidyr' was built under R version 4.1.2
#> Warning: package 'readr' was built under R version 4.1.2
#> Warning: package 'dplyr' was built under R version 4.1.2
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
```

We then take a look at the tutorialData inside PHONEMeS. It consists of
a top table results object from limma coming from phosphoproteomics.
This contains the information about phosphorylation sites that are
differntially abundant between two conditions.

``` r
head(tutorialData)
#>            ID     logFC  AveExpr         t      P.Value    adj.P.Val        B
#> 1    PSMA2_Y6  4.178268 16.69517  14.66493 1.719295e-30 3.374804e-26 58.57585
#> 2   HID1_S660 -3.051488 16.68615 -13.75250 9.323056e-28 7.803912e-24 52.37448
#> 3 HNRNPK_S116  5.490063 16.01733  13.61759 1.192712e-27 7.803912e-24 52.16460
#> 4    CAT_T511  3.283621 15.62845  13.51190 2.929472e-27 1.437565e-23 51.28095
#> 5 SH3BP5L_S30 -2.825266 16.82388 -13.29105 1.841500e-26 7.229361e-23 49.44288
#> 6   MDC1_T150  4.290834 15.96608  13.04989 2.819317e-26 9.223397e-23 49.06691
```

# Input preparation

This section provides a guide for the selection of deregulated
phosphorylation sites and kinases. It also includes an example of how to
estimate kinase activities using decoupleR which can be used for the
selection of deregulated kinases.

## Selection of deregulated phosphorylation sites

Deregulated phosphorylation sites can be selected based on their changes
detected within the differential analysis. We choose T values for the
selection over other metrices (p-values, logFC) as it combines the
significance and magnitude of the change. We then choose the top 10
deregulated phosphorylation sites that we will try to connect to
upstream kinases. This is in arbitrary number and users should decide on
a case to case basis what selection criteria to apply. The
phosphorylation sites could also be chosen based on a chosen percentage
or score cut-off.

``` r
top_pps <- tutorialData %>% 
  dplyr::filter(ID %in% phonemesPKN$target) %>%
  dplyr::arrange(dplyr::desc(base::abs(t))) %>%
  dplyr::slice(1:10)

deregulated_pps <- top_pps$t
names(deregulated_pps) <- top_pps$ID
```

## Selection of deregulated kinases

Kinases that are deregulated can either be selected manually, e.g. if a
specific kinase is known to be perturbed in an experiment, or based on
changes in their activity. Therefore, we will estimate kinase activity
changes using the decoupleR framework.

### Manual selection

For the manual selection, we create a names vector were kinases can be
classified as eiter up-regulated (1) or down-regulated (-1). In this
particular example we assume an up-regulation of AKT1.

``` r
deregulated_kinases_man <- c(AKT1 = 1)
```

### Kinase activity estimation

DecoupleR is a collection of computational methods used for
footprint-based analysis that assume omics data as a signature of
upstream biological activities. Here, we can use the changes in
phosphorylation sites as signatures for the activity of upstream
kinases. The interactions between phosphorylation sites and kinases are
again taken from the phonemesPKN. To use the phonemesPKN for decoupleR
we need to add a column for the likelihood of each interaction. Since
our network is not weighted, we just add a 1 for all interactions.

``` r
decoupler_network <- phonemesPKN %>% 
  dplyr::rename("mor" = interaction) %>% 
  tibble::add_column("likelihood" = 1)
```

As decoupleR functions expect a matrix as input, we are going to create
a one-column matrix from the tutorialData object. This one column matrix
contains the vector of T values obtained from limma. We again choose the
T values over other metrices.

``` r
decoupler_input <- tutorialData %>% 
  dplyr::filter(ID %in% decoupler_network$target) %>%
  tibble::column_to_rownames("ID") %>% 
  dplyr::select(t)
```

Before running any method in decoupleR, we need intersect our prior
knowledge network with the input matrix. With that we filter out
regulons with less than 5 target features.

``` r
decoupler_network <- decoupleR::intersect_regulons(mat = decoupler_input, 
                                                   network = decoupler_network, 
                                                   .source = source, 
                                                   .target = target, 
                                                   minsize = 5)
```

Additionally, the overlap of different regulons should be checked. We
filter out kinases with a high overlap in their target sets (correlation
\>= 0.9) and only keep one of them as these are affected by the multi
linear model used in mlm.

``` r
correlated_regulons <- decoupleR::check_corr(decoupler_network) %>% 
  dplyr::filter(correlation >= 0.9)

decoupler_network <- decoupler_network %>% 
  dplyr::filter(!source %in% correlated_regulons$source.2)
```

We use mlm to estimate kinase activity but if the user would like to
check other available methods they can run show_methods(). All methods
follow the same design pattern and arguments.

``` r
kinase_activity <- decoupleR::run_mlm(mat = decoupler_input, 
                                      network = decoupler_network)
```

We then choose the top 5 deregulated kinases that we will try to connect
to the downstream phosphorylation sites. This is again an arbitrary
number and users should decide on a case to case basis what selection
criteria to apply. Kinases could also be chosen based on a chosen
percentage or score cut-off. As before, kinases will eiter be classified
as up-regulated (1) or down-regulated (-1).

``` r
top_kinases <- kinase_activity %>% 
  dplyr::arrange(dplyr::desc(base::abs(score))) %>% 
  dplyr::slice(1:5)

deregulated_kinases <- top_kinases$score
names(deregulated_kinases) <- top_kinases$source

deregulated_kinases[deregulated_kinases > 0] <- 1
deregulated_kinases[deregulated_kinases < 0] <- -1
```

Moreover, we recommend the users to select kinases that do not appear to
be regulated based on their kinase activities. We assume that kinases
with a score below 0.5 are likely to not be deregulated between the
conditions and with that do not explain the observed changes. These
kinases will be filtered out from our network within PHONEMeS.

``` r
nc_kinases <- kinase_activity %>% 
  dplyr::filter(base::abs(score) <= 0.5) %>% 
  dplyr::pull(source)
```

# Run PHONEMeS

The deregulated phosphorylation sites and kinases can now be connected
using PHONEMeS, in specific the run_phonemes function. Here you can
either choose your manual kinase selection (deregulated_kinases_man) or
the kinases based on their activity changes (deregulated_kinases). For
simplicity reasons we will choose the manual kinase selection in this
example, as it will only take about 1 minute to be solved by PHONEMeS.
In the run_phonemes function the prior knowledge network will be pruned
by removing nodes 50 steps upstream and downstream of measurements and
inputs (default). Remember that the path to an interactive version of
IBM Cplex or CBC-COIN solver is required.

Before running PHONEMeS we need to set up the options for CARNIVAL which
is called within. The function default_carnival_options returns a list
with all variables in CARNIVAL based on your solver. For cplex and cbc
you will need to specify the solverPath.

``` r
carnival_options <- PHONEMeS::default_carnival_options(solver = "cplex")
carnival_options$solverPath <- "/Applications/CPLEX_Studio201/cplex/bin/x86-64_osx/cplex" #Example Path
```

``` r
phonemes_result <- PHONEMeS::run_phonemes(inputObj = deregulated_kinases_man, 
                                          measObj = deregulated_pps, 
                                          rmNodes = nc_kinases, 
                                          netObj = phonemesPKN,
                                          carnival_options = carnival_options)
#> Input nodes: 1 
#> Measurement nodes: 10 
#> Network nodes: 1537 
#> Network edges: 9289
#> --- Start of the CARNIVAL pipeline ---
#> 10:20:01 28.04.2022 Carnival flavour: vanilla
#> 10:20:01 28.04.2022 Generating variables for lp problem
#> 10:20:01 28.04.2022 Done: generating variables for lp problem
#> Saving preprocessed data.
#> Done: saving parsed data: /Users/smuellerdott/Documents/PHONEMeS/vignettes//parsedData_t10_20_01d28_04_2022n86.RData
#> 10:20:01 28.04.2022 Generating formulation for LP problem
#> 10:20:03 28.04.2022 Done: generating formulation for LP problem.
#> Saving LP file
#> Done: Saving LP file: /Users/smuellerdott/Documents/PHONEMeS/vignettes//lpFile_t10_20_01d28_04_2022n86.lp
#> 10:20:03 28.04.2022 Solving LP problem
#> Writing cplex command file
#> Done: writing cplex command file
#> Saving results...
#> 10:20:49 28.04.2022 Done: solving LP problem.
#> 10:20:49 28.04.2022 Getting the solution matrix
#> 10:20:50 28.04.2022 Done: getting the solution matrix.
#> 10:20:50 28.04.2022 Exporting solution matrix
#> 10:20:50 28.04.2022 Done: exporting solution matrix.
#> Cleaning intermediate files
#> Done: cleaning
#> 10:20:50 28.04.2022 All tasks finished.
#> 
#> --- End of the CARNIVAL pipeline ---
```

The output of PHONEMeS contains the network, measurements and inputs
used for the construction and the final results, including individual
solutions and one summary solution (weightedSIF, nodesAttributes). These
can be directly used to identify causal interactions between the
perturbed nodes and the selected kinases. In addition to extracting
direct information from the network, we can run different downstream
analysis based on the necessities of each project, e.g. Pathway
enrichment analysis.

## Re-attach phosphorylation sites

After the signaling network is constructed we can connect the
phosphorylation sites within the network to their protein if this is
also part of the network.

``` r
phonemes_result_pps <- PHONEMeS::reattach_psites(phonemes_result)
#> [1] "No psites to attach"
```

## Save network

The network sif and attributes file can be easily saved to a csv file.

``` r
readr::write_csv(phonemes_result_pps$res$weightedSIF, "examplePath/networkSIF.csv")
readr::write_csv(phonemes_result_pps$res$nodesAttributes, "examplePath/networkAtt.csv")
```

# Further functions

## Get protein network

Sometimes you might be interested in focusing on the protein
interactions of the network and removing phosphorylation sites. This can
also be helpful for the visualization of the network. You can extract
the protein network with the get_protein_network function. The output
will have the same structure as the phonemes_result object.

``` r
phonemes_result_protein <- get_protein_network(phonemes_result)
```

## Extract subnetwork

If you want to focus on a specific protein within the network you can
extract n-steps up- and/or down-stream of that node.

``` r
phonemes_result_sub <- PHONEMeS::extract_subnetwork(phonemes_result,
                                                    targets = "CDK2",
                                                    n = 1)
```

# Session Info

This tutorial was run on the date specified below.

``` r
Sys.Date()
#> [1] "2022-04-28"
```

The sessionInfo() at run time was:

``` r
sessionInfo()
#> R version 4.1.0 (2021-05-18)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Big Sur 10.16
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] forcats_0.5.1   stringr_1.4.0   dplyr_1.0.8     purrr_0.3.4    
#>  [5] readr_2.1.2     tidyr_1.2.0     tibble_3.1.6    ggplot2_3.3.5  
#>  [9] tidyverse_1.3.1 decoupleR_2.1.8 PHONEMeS_2.0.0 
#> 
#> loaded via a namespace (and not attached):
#>  [1] lubridate_1.8.0  lattice_0.20-45  assertthat_0.2.1 digest_0.6.29   
#>  [5] utf8_1.2.2       R6_2.5.1         cellranger_1.1.0 backports_1.4.1 
#>  [9] reprex_2.0.1     evaluate_0.15    httr_1.4.2       pillar_1.7.0    
#> [13] rlang_1.0.2      readxl_1.4.0     rstudioapi_0.13  Matrix_1.4-1    
#> [17] rmarkdown_2.14   bit_4.0.4        igraph_1.3.1     munsell_0.5.0   
#> [21] broom_0.8.0      compiler_4.1.0   modelr_0.1.8     xfun_0.30       
#> [25] pkgconfig_2.0.3  CARNIVAL_2.5.1   htmltools_0.5.2  tidyselect_1.1.2
#> [29] lpSolve_5.6.15   fansi_1.0.3      crayon_1.5.1     tzdb_0.3.0      
#> [33] dbplyr_2.1.1     withr_2.5.0      grid_4.1.0       jsonlite_1.8.0  
#> [37] gtable_0.3.0     lifecycle_1.0.1  DBI_1.1.2        magrittr_2.0.3  
#> [41] scales_1.2.0     cli_3.2.0        stringi_1.7.6    vroom_1.5.7     
#> [45] fs_1.5.2         xml2_1.3.3       ellipsis_0.3.2   generics_0.1.2  
#> [49] vctrs_0.4.1      rjson_0.2.21     tools_4.1.0      bit64_4.0.5     
#> [53] glue_1.6.2       hms_1.1.1        parallel_4.1.0   fastmap_1.1.0   
#> [57] yaml_2.3.5       colorspace_2.0-3 rvest_1.0.2      knitr_1.39      
#> [61] haven_2.5.0
```
