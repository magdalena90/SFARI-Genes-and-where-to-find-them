# SFARI Genes and where to find them; classification modelling to identify genes associated with Autism Spectrum Disorder from RNA-seq data

<br>

Code repository for **[SFARI Genes and where to find them; classification modelling to identify genes associated with Autism Spectrum Disorder from RNA-seq data](https://doi.org/10.1101/2021.01.29.428754)**

<br>

## Notes about this repository

All code is in R. The [drake](https://github.com/ropensci/drake) package is used to manage the workflow of the project, but the code can also be executed as a regular R script:

- The script `run.R` runs the project as a regular R script and saves the output in the **Results** folder

- The script `make.R` runs the project using `drake` and saves the output in the **.drake** folder, which can be accessed by name using `drake`'s `loadd()` function

<br>

**Note on using drake:** `Drake` provides a lot of useful features but it has two drawbacks in this project:

- Running the project using `run.R` is much faster than with `make.R` in computers with multiple cores because some packages use the `parallel` package underneath, which doesn't work well with `clustermq`, the package `drake` uses to distribute the work

- The Enrichment Analysis of the top modules is only available when running the code using `run.R`, because of compatibility issues between the package `clusterProfiler` and `drake`

<br>

## Running the code

1. Clone this repository

2. Download InputData from [doi.org/10.7488/ds/2980](https://doi.org/10.7488/ds/2980)

3. Execute `run.R` or `make.R` depending on whether you want your workflow to be run `drake` or not

<br>

## InputData

<br>

- **genes_GO_annotations:** Gene Ontology annotations for each gene

- **krishnan_probability_score.xlsx:** Krishnan's ASD probabilty score downloaded from [asd.princeton.edu](http://asd.princeton.edu/)

- **NCBI_gene2ensembl_20_02_07gz:** NCBI's mapping between genes symbols and ensembl IDs

- **NCBI_gene_info_20_02_07_.gz:** Functional annotations of the genes

- **RNAseq_ASD_datExpr.csv:** Gene expression matrix. Downloaded from [mgandal's github repository](https://github.com/mgandal/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap)

- **RNAseq_ASD_datMeta.csv:** Metadata of the samples from the gene expression matrix. Downloaded from [mgandal's github repository](https://github.com/mgandal/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap)

- **sanders_TADA_score.xlsx** Sanders TADA score downloaded from [He et al., 2013)[https://doi.org/10.1371/journal.pgen.1003671]

- **SFARI_genes_01-03-2020.csv:** SFARI Gene scores using new scoring system

- **SFARI_genes_08-29-2019.csv:** SFARI Gene scores using old scoring system

<br>

## Output

<br>

### Preprocessed Input Data

- **new_SFARI_dataset:** Dataframe with information about SFARI genes with the new annotation criteria (scores 1 to 3)

- **old_SFARI_dataset:** Dataframe with information about SFARI genes with the original annotation criteria (scores 1 to 6)

- **NCBI_dataset:** Dataframe with gene biotype annotation obtained from NCBI

- **GO_neuronal_dataset:** Dataframe with gene annotation indicating if they have some neuronal-related function in the Gene Ontology

- **Gandal_dataset:** RData object containing the preprocessed and normalised gene expression data

<br>

### WGCNA

- **modules_dataset:** Dataframe indicating the module each of the genes belong to

- **top_modules_by_Diagnosis:** Dataframe indicating the modules with the highest relation to Diagnosis as well as their correlation value

- **top_modules_by_SFARI**: Dataframe indicating the modules with the highest enrichment in SFARI Genes as well as their enrichment and adjusted p-value

- **top_modules_enrichment:** (not included in the `drake` workflow) Named list with the Enrichment results for all the modules with a strong correlation to Diagnosis or enriched in SFARI Genes

<br>

### Classification Model

- **classification_dataset:** Dataframe with the input data used for the classification models

- **biased_classification_model:** Named list with the information from the biased classification model, including the predictions for each gene and the coefficients and performance metrics of the model

- **unbiased_classification_model:** Named list with the information from the unbiased classification model. The list includes the same elements as biased_classification_model
