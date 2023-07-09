
# Microarray Data Analysis

This is an example of Microarray data analysis workflow.

This script uses **R programming language** to perform differential expression analysis on gene
expression data obtained from the NCBI Gene Expression Omnibus (GEO). The data used in this
script is from the GEO dataset **GSE10072** and contains gene expression values from 58 tumor
and 49 non-tumor tissues from 20 never smokers, 26 former smokers, and 28 current smokers.
The script loads several packages from the Bioconductor project and performs exploratory data
analysis (EDA) on the gene expression data before conducting differential expression analysis,
gene enrichment analysis and pathway analysis.


## Package installation:
Before using the script, several packages need to be installed, including BiocManager, limma,
sva, affyPLM, ggplot2, and preprocessCore.
The libraries required to be loaded are GEOquery, sva, limma, umap, maptools, clusterProfiler,
org.Hs.eg.db, ggplot2, dplyr


## Steps

- Data Loading
- Exploratory Data Analysis
- Preprocessing
- Surrogate Variable Analysis
- Normalization
- Differential Expression Analysis
- Gene Set Enrichment Analysis
- Pathway Analysis

For more information about individual steps refer to the Documentation.

## Data Visualization
Data is visualized at various stages using visualization tools such as

- Box plots
- Heatmaps
- Histograms
- Density plots
- QQ plots
- Volcano plots
- Bar plots

## Conclusion
Top 20 enriched pathways based on the gene expression data analyzed 
1. Vasculogenesis
2. Regulation of angiogenesis
3. Regulation of vasculature development
4. Positive regulation of angiogenesis
5. Positive regulation of vasculature development
6. Heart morphogenesis
7. Cell-substrate adhesion
8. Vascular process in circulatory system
9. Cardiac chamber morphogenesis
10. Regulation of cell-substrate adhesion
11. Regulation of focal adhesion assembly
12. Regulation of cell-substrate junction assembly
13. Endothelium development
14. Cell junction assembly
15. Cell-matrix adhesion
16. Regulation of cell-substrate junction organization
17. Smooth muscle cell differentiation
18. Regulation of cell-matrix adhesion
19. Tissue migration
20. Morphogenesis of a branching epithelium