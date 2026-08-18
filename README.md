# GSE148537-Differential-Expression-Analysis-of-MCF-7-shSPCA2-vs-Control (Microarray data)
This project performs differential gene expression analysis using the GSE148537 microarray dataset from the Gene Expression Omnibus (GEO).  
The dataset contains MCF-7 cells with:

shRNA control
shSPCA2 knockdown

The aim is to identify genes that are differentially expressed between shSPCA2 and control cells.

# Analysis Workflow

The analysis was performed in R using the following workflow:
Downloaded raw Affymetrix CEL files from GEO.
Read the CEL files using the affy package.
Performed RMA normalization.
Created experimental groups based on the sample phenotype.
Constructed a design matrix using limma.
Created the contrast:shSPCA2 - Control
Performed differential expression analysis using the limma linear model.
Applied empirical Bayes moderation using eBayes().
Annotated Affymetrix probe IDs using hgu133plus2.db.
Filtered significant DEGs using:
Adjusted P-value < 0.05
|log2 Fold Change| ≥ 1
Separated genes into upregulated and downregulated groups.
Saved the results as separate CSV files.

# Experimental Design

There were four biological samples:
MCF-7 shRNA control biological replicate 1	Control
MCF-7 shRNA control biological replicate 2	Control
MCF-7 shSPCA2 biological replicate 1	shSPCA2
MCF-7 shSPCA2 biological replicate 2	shSPCA2

# Differential Expression Results
shSPCA2VsControl = shSPCA2 - Control
Therefore:

Positive log2FC → higher expression in shSPCA2
Negative log2FC → lower expression in shSPCA2

Using the thresholds:

Adjusted P-value < 0.05
|log2FC| ≥ 1

# The analysis identified:

76 upregulated probe-level results
15 downregulated probe-level results
91 significant probe-level results in total

# R Packages
affy
GEOquery
limma
hgu133plus2.db
AnnotationDbi
tidyverse

# Next Steps

The next stages of the analysis will include:

Quality control and PCA
Sample correlation analysis
Heatmap of DEGs
Volcano plot
Gene Ontology (GO) enrichment
KEGG pathway enrichment
Biological interpretation of shSPCA2-associated transcriptional changes

Note: The experiment contains only two biological replicates per condition, so the results should be interpreted with consideration of the limited sample size.

