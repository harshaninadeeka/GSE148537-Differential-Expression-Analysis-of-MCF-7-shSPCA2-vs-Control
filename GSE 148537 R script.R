# script to perform RMA normalization
setwd("~/Training/Microarray")

library(affy)
library(GEOquery)
library(tidyverse)

# get supplementary files
getGEOSuppFiles("GSE148537")

# untar files
untar("GSE148537/GSE148537_RAW.tar", exdir = 'data1/')

# reading in .cel files
raw.data <- ReadAffy(celfile.path = "data1/")
raw.data

# performing RMA normalization
eset <- rma(raw.data)


# map probe IDs to gene symbols
gse <- getGEO("GSE148537", GSEMatrix = TRUE)
length(gse)

gse[[1]]@annotation
annotation(eset)
phenodata<-gse$GSE148537_series_matrix.txt.gz@phenoData
colnames(phenodata)

phenodata[, c("geo_accession", "title", "characteristics_ch1","status")]
phenodata$title



#Create the experimental groups
group <- factor(
  c("Control", "Control", "shSPCA2", "shSPCA2"),
  levels = c("Control", "shSPCA2")
)

group

#Install/load the annotation package

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hgu133plus2.db")


library(hgu133plus2.db)
library(AnnotationDbi)

#Load package
library(limma)
#Define design matrix
studyDesign = data.frame(Group = group)
designMatrix <- model.matrix(~ 0 + Group, data = studyDesign) 

colnames(designMatrix)[1:2] <- c("Control", "shSPCA2")

#Fit linear model
fit <- lmFit(exprs(eset), designMatrix)

#Define contrast matrix
contrastMatrix <- makeContrasts(shSPCA2VsControl = shSPCA2 - Control,
                                levels = designMatrix)
contrastMatrix

colnames(fit$coefficients)
ncol(fit$coefficients)
nrow(contrastMatrix)

#Specify contrasts
fit <- contrasts.fit(fit, contrastMatrix)




#Perform statistical testing
fit <- eBayes(fit, robust = TRUE)

#get all probes without filtering results
results <- topTable(fit, coef = "shSPCA2VsControl", number = nrow(eset),
                    p.value = 1)
print(dim(results))
print(head(results))


dim(exprs(eset))

#Convert probe IDs to gene symbols
probe_ids <- rownames(results)

gene_annotation <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = probe_ids,
  columns = c("SYMBOL", "ENTREZID", "GENENAME"),
  keytype = "PROBEID"
)
#check
head(gene_annotation)


#Join the gene IDs to limma results
results$PROBEID <- rownames(results)
results_annotated <- left_join(
  results,
  gene_annotation,
  by = "PROBEID"
)
head(results_annotated)

table(results_annotated$SYMBOL) |> 
  sort(decreasing = TRUE) |> 
  head(20)


#create DEG table (adj.P value <0.05)
DEGs <- results_annotated %>%
  filter(
    !is.na(SYMBOL),
    adj.P.Val < 0.05,
    abs(logFC) >= 1
  )
head(DEGs)
nrow(DEGs)




upregulated <- DEGs %>%
  filter(logFC >= 1)

downregulated <- DEGs %>%
  filter(logFC <= -1)
nrow(upregulated)
nrow(downregulated)



# Save complete differential expression results
write.csv(
  results_annotated,
  "GSE148537_all_DEG_results_annotated.csv",
  row.names = FALSE
)

# Save all significant DEGs
write.csv(
  DEGs,
  "GSE148537_DEGs_adjP0.05_logFC1.csv",
  row.names = FALSE
)

# Save upregulated genes
write.csv(
  upregulated,
  "GSE148537_upregulated_genes.csv",
  row.names = FALSE
)

# Save downregulated genes
write.csv(
  downregulated,
  "GSE148537_downregulated_genes.csv",
  row.names = FALSE
)

