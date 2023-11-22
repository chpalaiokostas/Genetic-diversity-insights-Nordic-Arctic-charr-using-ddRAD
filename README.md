# Genetic-diversity-insights-from-Nordic-Arctic-charr-using-ddRAD

The repository's code has been used to conduct a genetic diversity study on Nordic Arctic charr populations.

The study has been submitted to the scientific journal Aquaculture reports.

The following briefly describes the conducted data analysis:

## Variant calling

Processing of the sequenced reads, QC filtering and variant calling was performed through a Snakemake pipeline. The file `stacks2_env.yaml`contains detailed information about each software that snakemake used. Using this yaml file one should be able to replicate our analyses. Keep in mind that it is expected that snakemake and conda are already installed. Also the raw sequnece reads are publicly avalaible through NCBI under project.


## Custom population genetic analyses

The obtained vcf file was filtered using SNPfiltR (https://devonderaad.github.io/SNPfiltR/). Thereafter the R script `ac_scandinavia_pop_gen_analysis.R` was used for PCA, DAPC and for predicting population of origin. A confusion matrix plot can be made with Python code offerered in the Jupyter Notebook file `Confusion_matrix.ipynb`. 


## Detection of genetic clusters with machine learning

K-Means, Gaussian mixture and Bayesian Gaussian mixture models were fitted using Python's `scikit-learn` aiming to identify genetic clusters. The code Python code is available in the Jupyter Notebook file `ML_genetic_clusters/ML_genetic_clusters.ipynb`. 






