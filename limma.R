# Setup ------------------------------------------------------------------------

library("Biobase")
library("BiocGenerics")
library("limma")

rds <- "https://assets.datacamp.com/production/repositories/1626/datasets/0decef2850200efcf87b107b080959b31ec681ba/cll-eset.rds"
eset <- readRDS(url(rds))

x <- exprs(eset)
f <- fData(eset)
p <- pData(eset)

# Explore data -----------------------------------------------------------------

boxplot(x[1, ] ~ p[, "Disease"], main = f[1, "symbol"])
dim(eset)

# Subset to only include the first 10 samples (columns)
eset_sub <- eset[, 1:10]

# Check the dimensions of the subset
dim(eset_sub)

# Create a boxplot of the 1000th gene in eset_sub
boxplot(exprs(eset_sub)[1000, ] ~ pData(eset_sub)[, "Disease"], main = fData(eset_sub)[1000, "symbol"])

# limma pipeline ---------------------------------------------------------------

# Create design matrix for CLL study
design <- model.matrix(~Disease, data = pData(eset))

# Count the number of samples modeled by each coefficient
colSums(design)

# Load package
library(limma)

# Fit the model
fit <- lmFit(eset, design)

# Calculate the t-statistics
fit <- eBayes(fit)

# Summarize results
results <- decideTests(fit[, "Diseasestable"])
summary(results)

# group-means ------------------------------------------------------------------

# Create design matrix with no intercept
design <- model.matrix(~0 + Disease, data = pData(eset))

# Count the number of samples modeled by each coefficient
colSums(design)

# Load package
library(limma)

# Create a contrasts matrix
cm <- makeContrasts(status = Diseaseprogres. - Diseasestable,
                    levels = design)

# View the contrasts matrix
cm

# Load package
library(limma)

# Fit the model
fit <- lmFit(eset, design)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrasts = cm)

# Calculate the t-statistics for the contrasts
fit2 <- eBayes(fit2)

# Summarize results
results <- decideTests(fit2)
summary(results)

# Visualize results ------------------------------------------------------------

# Obtain the summary statistics for every gene
stats <- topTable(fit2, number = nrow(fit2), sort.by = "none")

# Plot a histogram of the p-values
hist(stats[, "P.Value"])

# Create a volcano plot. Highlight the top 5 genes
volcanoplot(fit2, highlight = 5, names = fit2$genes[, "symbol"])

# Enrichment -------------------------------------------------------------------

# Extract the entrez gene IDs
entrez <- fit2$genes[, "entrez"]

# Test for enriched GO categories
enrich_go <- goana(fit2, geneid = entrez, species = "Hs")

# View the top 20 enriched GO Biological Processes
topGO(enrich_go, ontology = "BP")

# Extract the entrez gene IDs
entrez <- fit2$genes[, "entrez"]

# Test for enriched KEGG Pathways
enrich_kegg <- kegga(fit2, geneid = entrez, species = "Hs")

# View the top 20 enriched KEGG pathways
topKEGG(enrich_kegg)
