# Introduction

Codes for generating the results in the paper: a statistical framework for differential pseudotime analysis in multiple single-cell RNA-seq 

by Wenpin Hou, Zhicheng Ji, Zeyu Chen, E. John Wherry, Stephanie C. Hicks\*, Hongkai Ji\*

## Folders 

**function**: main R functions used in this repository.

**h5func**: hdf5-related R functions.

**tree_variability**: module 1 analysis.

**hca_bone_marrow_data_analysis**: HCA-BM single-cell RNA-seq data analysis. The purpose of this HCA-BM case study is to illustrate all modules of Lamian and provide some benchmarking (2.2 - 2.4).

**hca_bone_marrow_data_integration_cutoff**: integration results using two other cutoffs for HCA-BM data.

**covid_data_analysis**: COVID19 data analysis. The purpose of this TB case study is to illustrate that Lamian can identify differential CD8 T cell transcriptional programs during a critical stage of disease severity transition (2.5).

**tb_data_analysis**: TB (tuberculosis) single-cell RNA-seq data analysis. The purpose of this TB  case study is to demonstrate and evaluate that Lamian can be applied to detect differences with respect to sample covariates while adjusting for batch effects (2.6). 

**harcohen_data_analysis**: Harcohen data analysis. For back up purpose only. This data is not presented.  

**age_differential**: identify age-associated genes.

**cellcycle_gene_analysis**: cellcycle gene analysis.

**efficiency**: computational efficiency.

## Citation 

A statistical framework for differential pseudotime analysis with multiple single-cell RNA-seq samples. 
Wenpin Hou, Zhicheng Ji, Zeyu Chen, E John Wherry, Stephanie C Hicks\*, Hongkai Ji\*. 
bioRxiv 2021.07.10.451910; doi: https://doi.org/10.1101/2021.07.10.451910. 

This manuscript is now under revision in a peer-review journal.

## Contact

Should you encounter any bugs or have any suggestions, please feel free to contact Wenpin Hou <wh2526@cumc.columbia.edu>, or open an issue on the Github page https://github.com/Winnie09/Lamian/issues.

