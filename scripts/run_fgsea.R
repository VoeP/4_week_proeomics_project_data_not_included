
### Following example from 
# https://davetang.org/muse/2018/01/10/using-fast-preranked-gene-set-enrichment-analysis-fgsea-package/


# R package fgsea
library(fgsea)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("reactome.db")

# Load a list of significant genes from the study GSE14308 (Mus musculus)
# Wei et al, Immunity. 2009 January 16; 30(1): 155–167. doi:10.1016/j.immuni.2008.12.009.

# Global Mapping of H3K4me3 and H3K27me3 Reveals Specificity
# and Plasticity in Lineage Fate Determination of Differentiating
# CD4 + T Cells

ranks <- read.table("rnk_GEO_GSE14308.file",
                  header=TRUE, colClasses = c("character", "numeric"))
ranks <- setNames(ranks$t, ranks$rn)


# We will check gene set enrichment at the top of the rank list 
# with gene sets from the reactome pathways database

library(reactome.db)
reactome_pathways <- reactomePathways(names(ranks))


# Run fgsea:
fgsea_reactome <- fgsea(pathways = reactome_pathways, 
                        stats = ranks,
                        minSize=15,
                        maxSize=500,
                        nperm=100000)

# order fgsea_reactome in order of decreasing pathways signicance
fgsea_reactome <-fgsea_reactome[order(pval),]

# Top 10 pathways
head(fgsea_reactome,10)

# All significant pathways
significant_pathways<-fgsea_reactome[fgsea_reactome$padj<=0.01,]



# Visualise enrichment score for top pathway
plotEnrichment(reactome_pathways[["Cell Cycle"]],ranks) 

# Visualise enrichment score for last significant pathway
plotEnrichment(reactome_pathways[["Kinesins"]],ranks)

# Visualise enrichment score for not significant pathways:

plotEnrichment(reactome_pathways[["GTP hydrolysis and joining of the 60S ribosomal subunit"]],ranks)

plotEnrichment(reactome_pathways[["Complex I biogenesis"]],ranks)
plotEnrichment(reactome_pathways[["TBC/RABGAPs"]],ranks)
plotEnrichment(reactome_pathways[["RET signaling"]],ranks)



######Construction from duplicate####




####construction
mean()


# R package fgsea
library(fgsea)

#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("reactome.db")

# Load a list of significant genes from the study GSE14308 (Mus musculus)
# Wei et al, Immunity. 2009 January 16; 30(1): 155–167. doi:10.1016/j.immuni.2008.12.009.

# Global Mapping of H3K4me3 and H3K27me3 Reveals Specificity
# and Plasticity in Lineage Fate Determination of Differentiating
# CD4 + T Cells

ranks <- read.table("rnk_GEO_GSE14308.file",
                    header=TRUE, colClasses = c("character", "numeric"))
ranks <- setNames(ranks$t, ranks$rn)


# We will check gene set enrichment at the top of the rank list 
# with gene sets from the reactome pathways database

library(reactome.db)
reactome_pathways <- reactomePathways(names(ranks))


# Run fgsea:
fgsea_reactome <- fgsea(pathways = reactome_pathways, 
                        stats = ranks,
                        minSize=15,
                        maxSize=500,
                        nperm=100000)

# order fgsea_reactome in order of decreasing pathways signicance
fgsea_reactome <-fgsea_reactome[order(pval),]

# Top 10 pathways
head(fgsea_reactome,10)

# All significant pathways
significant_pathways<-fgsea_reactome[fgsea_reactome$padj<=0.01,]



# Visualise enrichment score for top pathway
plotEnrichment(reactome_pathways[["Cell Cycle"]],ranks) 

# Visualise enrichment score for last significant pathway
plotEnrichment(reactome_pathways[["Kinesins"]],ranks)

# Visualise enrichment score for not significant pathways:

plotEnrichment(reactome_pathways[["GTP hydrolysis and joining of the 60S ribosomal subunit"]],ranks)

plotEnrichment(reactome_pathways[["Complex I biogenesis"]],ranks)
plotEnrichment(reactome_pathways[["TBC/RABGAPs"]],ranks)
plotEnrichment(reactome_pathways[["RET signaling"]],ranks)
