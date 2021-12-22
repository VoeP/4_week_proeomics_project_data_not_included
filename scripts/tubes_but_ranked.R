####ranked analysis of tubes

library(openxlsx)

PFP_TS<-read.xlsx("EV_352_PFP_TS-C.xlsx")
PFP_C<-read.xlsx("EV_PFP_C_anova.xlsx")




library(SetRankCorr)
library(dplyr)
library(tidyr)
library(Homo.sapiens)
library(EnrichmentBrowser)
library(SetRankCorr)


uniprot2EntrezID=  createIDConverter("Homo.sapiens", "UNIPROT", "ENTREZID")


ExtractACs <- function(ProtACs){
  mylistsemicol=unlist(strsplit(as.character(ProtACs), ";"))
  mylist=sapply(strsplit(as.character(mylistsemicol), "\\|"), "[[", 2)
  # take canonical form:
  canonicals<-sapply(mylist,function(x) strsplit(x,"-")[[1]][1])
  # eliminate duplicates:
  canonicals<-unique(canonicals)
  myoutputstring=paste(canonicals,collapse = ",") 
  return(myoutputstring)
}


ProtENTREZ_bygroup <- function(ProtACs,uniprot2Entrez){
  cat(ProtACs,"\n") # detect possible problems
  mylistsemicol=unlist(strsplit(as.character(ProtACs), ",")) 
  myENTREZlist=uniprot2Entrez(na.omit(mylistsemicol),
                              na.rm=TRUE,drop.ambiguous = TRUE)
  myoutputstring=NA
  myoutputstring=paste(na.omit(unique(myENTREZlist)),collapse = ",")
  return(myoutputstring)
}


PFP_TS$ProtsACs=sapply(PFP_TS[,"Protein.IDs"],ExtractACs)
PFP_TS<-PFP_TS[!is.na(PFP_TS[,"iiTop3.tt_adjpVal.352_PFP_TS.-.352_PFP_C"]),]

PFP_TS=PFP_TS %>% 
  mutate(ProtsACs= strsplit(as.character(ProtsACs), ",")) %>% 
  unnest(ProtsACs)
PFP_TS$ENTREZ<-sapply(PFP_TS$ProtsACs,ProtENTREZ_bygroup,uniprot2Entrez=uniprot2EntrezID)
PFP_TS_noENTREZ<-PFP_TS[which(PFP_TS$ENTREZ==""),]
PFP_TS<-PFP_TS[which(PFP_TS$ENTREZ != ""),]





PFP_C$ProtsACs=sapply(PFP_C[,"Protein.group"],ExtractACs)
PFP_C<-PFP_C[!is.na(PFP_C[,"log.ANOVA.p.value"]),]

PFP_C=PFP_C %>% 
  mutate(Protein.group= strsplit(as.character(Protein.group), ";")) %>% 
  unnest(Protein.group)
PFP_C$ENTREZ<-sapply(PFP_C$Protein.group,ProtENTREZ_bygroup,uniprot2Entrez=uniprot2EntrezID)
PFP_C_noENTREZ<-PFP_C[which(PFP_C$ENTREZ==""),]
PFP_C<-PFP_C[which(PFP_C$ENTREZ != ""),]



library(GeneSets.Homo.sapiens)

options(mc.cores=2)
TS_reactom<- buildSetCollection(REACTOME, referenceSet=PFP_TS$ENTREZ, maxSetSize = 500)
PFPC_reactome<- buildSetCollection(REACTOME, referenceSet=PFP_C$ENTREZ, maxSetSize = 500)

save(TS_reactom, file="TS_reactome.RData")
save(PFPC_reactome, file="PFPC_reactome.RData")


options(mc.cores=2)
library(SetRankCorr)

network_ranked_tube = setRankAnalysis(PFP_TS$ENTREZ, TS_reactom, use.ranks=TRUE, setPCutoff = 0.01,
                                      fdrCutoff = 0.05)

exportSingleResult(network_ranked_tube, PFP_TS$ENTREZ, TS_reactom, networkName = "tube_system_network_ranked", IDConverter=NULL)



network_ranked_tube = setRankAnalysis(PFP_C$ENTREZ, PFPC_reactome, use.ranks=TRUE, setPCutoff = 0.01,
                                      fdrCutoff = 0.05)

exportSingleResult(network_ranked_tube, PFP_C$ENTREZ, PFPC_reactome, networkName = "carrier_system_network_ranked", IDConverter=NULL)






####fgsea for tube system


attach(PFP_TS)
PFP_TS$fcSign=sign(`iiTop3.log2FC.352_PPP_TS.-.352_PFP_TS`)

PFP_TS$metric=-(log10(PFP_TS$`iiTop3.tt_adjpVal.352_PFP_TS.-.352_PFP_C`)) / PFP_TS$fcSign
PFP_TS<-PFP_TS[order(PFP_TS$metric),]
y<-PFP_TS[,c("ENTREZ", "metric")]
write.table(y,file="PFP_TS_expression.rnk",quote=F,sep="\t",row.names=F)



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

ranks <- read.table("PFP_TS_expression.rnk",
                    header=TRUE, colClasses = c("character", "numeric"))
ranks <- setNames(ranks$metric, ranks$ENTREZ)


# We will check gene set enrichment at the top of the rank list 
# with gene sets from the reactome pathways database

library(reactome.db)
reactome_pathways <- reactomePathways(names(ranks))


# Run fgsea:
fgsea_PFP_TS <- fgsea(pathways = reactome_pathways, 
                   stats = ranks,
                   minSize=15,
                   maxSize=500
)

# order fgsea_reactome in order of decreasing pathways signicance
fgsea_PFP_TS <-fgsea_PFP_TS[order(pval),]

# Top 10 pathways
head(fgsea_PFP_TS,10)

# All significant pathways
significant_pathways_PFP_TS<-fgsea_PFP_TS[fgsea_PFP_TS$padj<=0.01,]

plotEnrichment(reactome_pathways[["B-WICH complex positively regulates rRNA expression"]],ranks)






####fgsea for carrier

attach(PFP_C)
## fold change doesn't exist; PFP_C$fcSign=sign(`iiTop3.log2FC.352_PPP_TS.-.352_PFP_TS`)

PFP_C$metric=-PFP_C$`-Log.ANOVA.p.value`
PFP_C<-PFP_C[order(PFP_C$metric),]
y<-PFP_C[,c("ENTREZ", "metric")]
write.table(y,file="PFP_C_expression.rnk",quote=F,sep="\t",row.names=F)



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

ranks <- read.table("PFP_C_expression.rnk",
                    header=TRUE, colClasses = c("character", "numeric"))
ranks <- setNames(ranks$metric, ranks$ENTREZ)


# We will check gene set enrichment at the top of the rank list 
# with gene sets from the reactome pathways database

library(reactome.db)
reactome_pathways <- reactomePathways(names(ranks))


# Run fgsea:
fgsea_PFP_C <- fgsea(pathways = reactome_pathways, 
                      stats = ranks,
                      minSize=15,
                      maxSize=500
)

# order fgsea_reactome in order of decreasing pathways signicance
fgsea_PFP_C <-fgsea_PFP_C[order(pval),]

# Top 10 pathways
head(fgsea_PFP_C,10)

# All significant pathways
significant_pathways_PFP_C<-fgsea_PFP_C[fgsea_PFP_C$padj<=0.01,]





library(gprofiler2)
gostrePFP_TS<-gost(
  query=PFP_TS$ENTREZ,
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

gostrePFP_C<-gost(
  query=PFP_C$ENTREZ,
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

panther_pfpc<-data.frame(I(PFP_C$ENTREZ), I(PFP_C$metric))
write.table(panther_pfpc, file="pfp_c_panther_input.txt", quote=FALSE, row.names = FALSE, sep = "\t")

panther_pfpts<-data.frame(I(PFP_TS$ENTREZ), I(PFP_TS$metric))
write.table(panther_pfpts, file="pfp_ts_panther_input.txt", quote=FALSE, row.names = FALSE, sep="\t")



gost_pfpts<-data.frame(gostrePFP_TS[["result"]][["term_name"]], gostrePFP_TS[["result"]][["source"]])
gost_combined_PFPTS_REAC<-gost_pfpts[which(gost_pfpts[[2]]=="REAC"),]


fgsea_names_pfpts<- fgsea_PFP_TS$pathway
gost_ts<-gost_combined_PFPTS_REAC$gostrePFP_TS...result......term_name...
panther_ts_path<-read.table("panther_output_pfp_ts.txt", sep = "\t", col.names = paste0("V",seq_len(10)), fill = TRUE)
panther_ts_path<-panther_ts_path[-1,]
panther_ts_path$V1<-gsub("\\(R.*)$","",panther_ts_path$V1)
panther_ts_path$V1<-trimws(panther_ts_path$V1, which="right")

colors <- c("#6b7fff", "#c3db0f", "#FF0000", "#FF6666")

intersect(gost_ts, fgsea_names_pfpts)


library(ggVennDiagram)
library(VennDiagram)
venn.diagram(x = list(fgsea_names_pfpts, gost_ts, panther_ts_path$V1, setrank_pfpts$V2) ,
             category.names = c("fgsea", "gost", "panther", "SetRank"),
             filename = 'pfpts_fgsea_with_setrankvenn.png',
             output=TRUE,
             imagetype="png", 
             scaled = FALSE,
             col = "black",
             fill = colors,
             cat.col = colors,
             cat.cex = 2,
             margin = 0.15
)


first_intersect<-intersect(gost_ts, fgsea_names_pfpts)
second_intersect<-intersect(first_intersect, panther_ts_path$V1)
intersect(significant_pathways_PFP_TS, second_intersect)
intersect(first_intersect, significant_pathways_PFP_TS$pathway)
setrank_pfpts<-read.table("tube_system_network_ranked_pathways.txt", sep="\t")
intersect(second_intersect, setrank_pfpts$V2)

setr_c<-read.table("carrier_system_network_ranked_pathways.txt", fill = TRUE, sep = "\t", header = TRUE)

