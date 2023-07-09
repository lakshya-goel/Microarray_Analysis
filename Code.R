if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install('limma',force = TRUE)
BiocManager::install('sva',force = TRUE)
BiocManager::install('affyPLM',force = TRUE)
BiocManager::install("ggplot2", dependencies = TRUE,force=TRUE)
BiocManager::install('biomaRt',force = TRUE)

library(GEOquery)
library(sva)
library(limma)
library(umap)
library(maptools)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

gse <- getGEO("GSE10072", GSEMatrix = TRUE,AnnotGPL=TRUE) #107 final expression values from 58 tumor and 49 non-tumor tissues from 20 never smokers, 26 former smokers, and 28 current smokers.
gse <- gse[[1]]
exprs <- exprs(gse)  #matrix of gene expression values with rows representing genes and columns representing samples.

#EDA
dim(exprs)                 #dataset has 22283 features and 107 samples
boxplot(exprs, outline=TRUE)  #This will display boxplots of gene expression values with outliers indicated. Outlier samples can be due to technical issues or biological factors and may need to be removed or treated separately in downstream analysis.
heatmap(cor(exprs))    #heatmap of correlation between samples
sum(is.na(exprs))    #This will display the total number of missing values in the gene expression data. Missing values can affect downstream analysis and may need to be imputed or removed.

pdata <- pData(gse)  #a data frame with information about the samples, including their class labels.
design <- model.matrix(~pData(gse)$characteristics_ch1.1)
mod <- sva(exprs,design)            #Number of significant surrogate variables is:  11 
                                    #Iteration (out of 5 ):1  2  3  4  5  

means <- apply(exprs, 1, mean)
vars <- apply(exprs, 1, var)

# Load necessary packages
exprs <- normalizeBetweenArrays(exprs(gse))
pdata <- pData(gse)
fdata <- fData(gse)

# List the data attributes
attributes(gse)

##log transforming the expression values and plotting density plot and histogram to see the difference

# Log2 transformation of expression values
exprs_log2 <- log2(exprs)

# Plot density plot of original and log-transformed data
par(mfrow = c(1, 2))
plot(density(exprs), main = "Density Plot of Original Data")
plot(density(exprs_log2), main = "Density Plot of Log-Transformed Data")

# Plot histogram of original data
hist(exprs, main = "Histogram of Original Data")
# Plot histogram of log-transformed data
hist(exprs_log2, main = "Histogram of Log-Transformed Data")

# group membership for all samples
gsms <- paste0("01010001011001010011010010101010010101001010101010",
               "01101010101000101010100010011001010111101010100010",
               "1010011")
sml <- strsplit(gsms, split="")[[1]]

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Lung_Cancer","Healthy_Control"))
levels(gs) <- groups
gse$group <- gs
design <- model.matrix(~group + 0, gse)
colnames(design) <- levels(gs)

ord <- order(gs) # order samples by group
palette(c(
    "#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
    "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"
))
par(mfrow = c(1, 2))
#par(mar = c(7, 4, 2, 1))
title1 <- paste("GSE10072", "/", annotation(gse), sep = "")
boxplot(exprs[, ord], boxwex = 0.6, notch = F, main = title1, outline = FALSE, las = 2, col = gs[ord])
legend("topleft", groups, fill = palette(), bty = "n")

#par(mar = c(7, 4, 2, 1))
title2 <- paste("GSE10072 (log transformed)", "/", annotation(gse), sep = "")
boxplot(exprs_log2[, ord], boxwex = 0.6, notch = F, main = title2, outline = FALSE, las = 2, col = gs[ord])
legend("topleft", groups, fill = palette(), bty = "n")
dev.off()

#differential expression analysis
exprs(gse) <- exprs_log2

# make proper column names to match toptable 
fLabels <- make.names(fvarLabels(gse))

#t test with holm's correction
healthy <- which(pData(gse)$group == "Healthy_Control") 
cancer <- which(pData(gse)$group == "Lung_Cancer")
ttest <- apply(exprs, 1, function(x) t.test(x[healthy], x[cancer]))
logfc <- apply(exprs, 1, function(x) mean(x[cancer]) - mean(x[healthy]))
pval <- sapply(ttest, function(x) x$p.value)
adj_pval <- p.adjust(pval, method = "holm")
results <- data.frame(logFC = logfc, PValue = pval, Adj_PValue = adj_pval)

# volcano plot
ggplot(results, aes(x = logFC, y = -log10(PValue))) +
    geom_point(aes(color = Adj_PValue < 0.05 & abs(logfc)<1), size = 1) +
    scale_color_manual(values = c("blue", "red")) +
    theme_classic() +
    ggtitle("Volcano Plot") +
    xlab("Log2(FC)") +
    ylab("-Log10(p-value)")

fit <- lmFit(exprs, design)  # fit linear model
# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="holm", sort.by="B", number=250)  #The table includes up to 250 genes sorted by the absolute value of the log-odds ratio (B).
gene_list <- rownames(tT)
tT$ID <- gene_list
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC"))
# Visualize and quality control test results.
# Build histogram of P-values for all genes.

tT2 <- topTable(fit2, adjust="holm", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="holm", p.value=0.05)


# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) #identifying "good" probes for the Q-Q plot. 
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gse)

# box-and-whisker plot
dev.new(width=3+ncol(gse)/6, height=5)
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE10072", "/", annotation(gse), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
dev.off()

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE10072", "/", annotation(gse), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)

pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE10072")

expr_mat <- exprs(gse)

# Load required package
library(biomaRt)

# Specify the mart to use (e.g. for human)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="useast.ensembl.org")


# Get the gene symbols from probe IDs
gene_list <- getBM(
    attributes = c("entrezgene_id", "affy_hg_u133_plus_2"),
    filters = "affy_hg_u133_plus_2",
    values = gene_list,
    mart = mart
)$entrezgene_id

# Run GO enrichment analysis with gene symbols
ego <- enrichGO(gene = gene_list,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.1,
                readable = TRUE)

# Show the enriched GO terms
ego@result

#Pathways
pathways <- head(ego@result, 20)
pathways$Description

#visualization
dot_plot <- ggplot(pathways, aes(x = -log10(pvalue), y = Description))
dot_plot + geom_point(size = 2) + xlab("-log10(pvalue)") + ylab("Pathway") 

bar_plot <- ggplot(pathways, aes(y = reorder(Description, -log10(pvalue)), x = -log10(pvalue))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    ylab("Pathway") +
    xlab("-log10(p-value)") +
    ggtitle("Top 20 enriched pathways")
bar_plot + theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
