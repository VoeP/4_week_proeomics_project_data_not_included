library(GSEABenchmarkeR)
library(SetRankCorr)
library(EnrichmentBrowser)
library(Homo.sapiens)
library(openxlsx)
library(gprofiler2)


geo2kegg <- loadEData("geo2kegg")

go.gs <- getGenesets(org="hsa", db="go")
go.ora.res <- runEA(geo2kegg[27:28], method="ora", gs=go.gs, perm=0)




geo2kegg <- loadEData("geo2kegg")
geo2kegg <- maPreproc(geo2kegg)
geo2kegg <- runDE(geo2kegg, de.method="limma", padj.method="flexible")
go.gs <- getGenesets(org="hsa", db="go")
go.ora.res <- runEA(geo2kegg[[27]], method=c("padog"), gs=go.gs,
                    perm=1000)

pado<-as.data.frame(go.ora.res$padog)
pado<-pado[order(pado$GSE14924_CD4.ranking.PVAL),]
pado_trim<-pado[which(pado$GSE14924_CD4.ranking.PVAL<0.05),]
write(pado_trim$GSE14924_CD4.ranking.GENE.SET, "pado_otput_CD4.txt", ncolumns = 1)




go.ora.res2 <- runEA(geo2kegg[[27]], method=c("ora"), gs=go.gs,
                     perm=0)

oraRes<-as.data.frame(go.ora.res2$ora)
oraRes<- oraRes[order(oraRes$GSE14924_CD4.ranking.PVAL)]
ora_trim<-oraRes[which(oraRes$GSE14924_CD4.ranking.PVAL<0.05),]
write(ora_trim$GSE14924_CD4.ranking.GENE.SET, "ora_otput_CD4.txt", ncolumns = 1)






go.ora.res3 <- runEA(geo2kegg[[27]], method=c("camera"), gs=go.gs,
                     perm=0)


cam<-as.data.frame(go.ora.res3$camera)
cam<-cam[order(cam$GSE14924_CD4.ranking.PVAL)]
cam_trim<-cam[which(cam$GSE14924_CD4.ranking.PVAL<0.05),]
write(cam_trim$GSE14924_CD4.ranking.GENE.SET, "cam_otput_CD4.txt", ncolumns = 1)




kegg.gs <- getGenesets(org="hsa", db="kegg")
kegg.ora.res <- runEA(geo2kegg[[1]], method="ora", gs=kegg.gs, perm=0)
kegg.ora.res




getwd()

cd4_dataset<-read.table("GSE14924_CD4_ENTREZ_ordered.csv", header=TRUE)
cd4_dataset<-cd4_dataset[order(cd4_dataset$padj),]
cd4_dataset_foreground<-subset(cd4_dataset, padj<0.05)

write(cd4_dataset_foreground$ENTREZID, "cd4_fore.txt", ncolumns = 1)

##honestly, the data might just be too large.-->ask what to do if the data is too large.
gost_cd4<-gost(
  query=unique(na.omit(cd4_dataset$ENTREZID)),
  organism = "hsapiens",
  ordered_query = TRUE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "ENTREZGENE_ACC",
  sources = NULL,
  as_short_link = FALSE
)






library(ggVennDiagram)
library(VennDiagram)

cam_CD4<-read.table("cam_CD4.txt")
ora_CD4<-read.table("ora_CD4.txt")
pado_CD4<-read.table("pado_CD4.txt")

##this line just proves that all identifiers were unique: GOBP_uniq<-unique(GOBP$V1), since same nr of entries


colors <- c("#6b7fff", "#c3db0f", "#FF0000")


venn.diagram(x = list(cam_CD4$V1, ora_CD4$V1, pado_CD4$V1) ,
             category.names = c("camera", "ora", "padog"),
             filename = 'cd4_padooracam_venn.png',
             output=TRUE,
             imagetype="png", 
             scaled = FALSE,
             col = "black",
             fill = colors,
             cat.col = colors,
             cat.cex = 2,
             margin = 0.15
)




gostreCD4<-gost(
  query=cd4_dataset$ENTREZID,
  organism = "hsapiens",
  ordered_query = TRUE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c( "custom"),
  custom_bg = NULL,
  numeric_ns = "ENTREZGENE_ACC",
  sources = NULL,
  as_short_link = FALSE
)

library(openxlsx)



#######FGSEA for CD4 


attach(cd4_dataset)
cd4_dataset$fcSign=sign(`log2FC_d_c`)

cd4_dataset$metric=cd4_dataset$padj / cd4_dataset$fcSign
cd4_dataset<-cd4_dataset[order(cd4_dataset$metric),]
y<-cd4_dataset[,c("ENTREZID", "metric")]
write.table(y,file="expression_cd4.rnk",quote=F,sep="\t",row.names=F)


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

ranks <- read.table("expression_cd4.rnk",
                    header=TRUE, colClasses = c("character", "numeric"))
ranks <- setNames(ranks$metric, ranks$ENTREZID)


# We will check gene set enrichment at the top of the rank list 
# with gene sets from the reactome pathways database

library(reactome.db)
reactome_pathways <- reactomePathways(names(ranks))


# Run fgsea:
fgsea_cd4 <- fgsea(pathways = reactome_pathways, 
                   stats = ranks,
                   minSize=15,
                   maxSize=500
)

# order fgsea_reactome in order of decreasing pathways signicance
fgsea_cd4 <-fgsea_cd4[order(pval),]

# Top 10 pathways
head(fgsea_cd4,10)

# All significant pathways
significant_pathways_cd4<-fgsea_cd4[fgsea_cd4$padj<=0.01,]








####attempt at getting reac stuff instead of GO



geo2kegg <- loadEData("geo2kegg")

go.gs <- getGenesets(org="hsa", db="reac")
go.ora.res <- runEA(geo2kegg[27:28], method="ora", gs=go.gs, perm=0)




geo2kegg <- loadEData("geo2kegg")
geo2kegg <- maPreproc(geo2kegg)
geo2kegg <- runDE(geo2kegg, de.method="limma", padj.method="flexible")
go.gs <- getGenesets(org="hsa", db="go")
go.ora.res <- runEA(geo2kegg[[27]], method=c("padog"), gs=go.gs,
                    perm=1000)

pado<-as.data.frame(go.ora.res$padog)
pado<-pado[order(pado$GSE14924_CD4.ranking.PVAL),]
pado_trim<-pado[which(pado$GSE14924_CD4.ranking.PVAL<0.05),]
write(pado_trim$GSE14924_CD4.ranking.GENE.SET, "pado_otput_CD4.txt", ncolumns = 1)




go.ora.res2 <- runEA(geo2kegg[[27]], method=c("ora"), gs=go.gs,
                     perm=0)

oraRes<-as.data.frame(go.ora.res2$ora)
oraRes<- oraRes[order(oraRes$GSE14924_CD4.ranking.PVAL)]
ora_trim<-oraRes[which(oraRes$GSE14924_CD4.ranking.PVAL<0.05),]
write(ora_trim$GSE14924_CD4.ranking.GENE.SET, "ora_otput_CD4.txt", ncolumns = 1)






go.ora.res3 <- runEA(geo2kegg[[27]], method=c("camera"), gs=go.gs,
                     perm=0)


cam<-as.data.frame(go.ora.res3$camera)
cam<-cam[order(cam$GSE14924_CD4.ranking.PVAL)]
cam_trim<-cam[which(cam$GSE14924_CD4.ranking.PVAL<0.05),]
write(cam_trim$GSE14924_CD4.ranking.GENE.SET, "cam_otput_CD4.txt", ncolumns = 1)




kegg.gs <- getGenesets(org="hsa", db="kegg")
kegg.ora.res <- runEA(geo2kegg[[1]], method="ora", gs=kegg.gs, perm=0)
kegg.ora.res




getwd()

cd4_dataset<-read.table("GSE14924_CD4_ENTREZ_ordered.csv", header=TRUE)
cd4_dataset<-cd4_dataset[order(cd4_dataset$padj),]
cd4_dataset_foreground<-subset(cd4_dataset, padj<0.05)

write(cd4_dataset_foreground$ENTREZID, "cd4_fore.txt", ncolumns = 1)

##honestly, the data might just be too large.-->ask what to do if the data is too large.
gost_cd4<-gost(
  query=unique(na.omit(cd4_dataset$ENTREZID)),
  organism = "hsapiens",
  ordered_query = TRUE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS",
                        "analytical"),
  domain_scope = c("annotated", "known", "custom", "custom_annotated"),
  custom_bg = NULL,
  numeric_ns = "ENTREZGENE_ACC",
  sources = NULL,
  as_short_link = FALSE
)




