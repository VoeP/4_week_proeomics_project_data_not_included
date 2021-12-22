####Huh7 ranked, I guess?

library(openxlsx)

Huh7_cellular<-read.xlsx("Huh-7 cellular proteome.xlsx")
Huh7_cellular_foreground<-subset(Huh7_cellular, `DV2_Con.-Log10.Student's.T-test.p-value`>1.3)


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



Huh7_cellular$ENTREZ<-sapply(Huh7_cellular$Accession,ProtENTREZ_bygroup,uniprot2Entrez=uniprot2EntrezID)
Huh7_noENTREZ<-Huh7_cellular[which(Huh7_cellular$ENTREZ==""),]
Huh7_cellular<-Huh7_cellular[which(Huh7_cellular$ENTREZ != ""),]

Huh7_cellular_foreground$ENTREZ<-sapply(Huh7_cellular_foreground$Accession,ProtENTREZ_bygroup,uniprot2Entrez=uniprot2EntrezID)
Huh7_foreground_noENTREZ<-Huh7_cellular_foreground[which(Huh7_cellular_foreground$ENTREZ==""),]
Huh7_cellular_foreground<-Huh7_cellular_foreground[which(Huh7_cellular_foreground$ENTREZ != ""),]


write(Huh7_cellular$ENTREZ, "background_Huh7_cellular.txt", ncolumns = 1)
write(Huh7_cellular_foreground$ENTREZ, "foreground_Huh7_cellular.txt", ncolumns = 1)



Huh7_secretome<-read.xlsx("Huh-7 cellular secretome.xlsx")
Huh7_secretome_foreground<-subset(Huh7_cellular, `DV2_Con.-Log10.Student's.T-test.p-value`>1.3)


Huh7_secretome$ENTREZ<-sapply(Huh7_secretome$Accession,ProtENTREZ_bygroup,uniprot2Entrez=uniprot2EntrezID)
Huh7_secretome_noENTREZ<-Huh7_secretome[which(Huh7_secretome$ENTREZ==""),]
Huh7_secretome<-Huh7_secretome[which(Huh7_secretome$ENTREZ != ""),]

Huh7_secretome_foreground$ENTREZ<-sapply(Huh7_secretome_foreground$Accession,ProtENTREZ_bygroup,uniprot2Entrez=uniprot2EntrezID)
Huh7_secretome__foreground_noENTREZ<-Huh7_secretome_foreground[which(Huh7_secretome_foreground$ENTREZ==""),]
Huh7_secretome_foreground<-Huh7_secretome_foreground[which(Huh7_secretome_foreground$ENTREZ != ""),]

write(Huh7_secretome$ENTREZ, "background_Huh7_secretome.txt", ncolumns = 1)
write(Huh7_secretome_foreground$ENTREZ, "foreground_Huh7_secretome.txt", ncolumns = 1)

library(GeneSets.Homo.sapiens)
library(GSEABenchmarkeR)
library(EnrichmentBrowser)
library(annotate)
library(SetRankCorr)

load("collection_Huh7_proteome_500_GOBP.RData")
load("collection_Huh7_secretome_500_GOBP.RData")


Huh7_cellular_foreground<-Huh7_cellular_foreground[order(-(Huh7_cellular_foreground$`DV2_Con.-Log10.Student's.T-test.p-value`)),]
Huh7_secretome_foreground<-Huh7_secretome_foreground[order(-(Huh7_cellular_foreground$`DV2_Con.-Log10.Student's.T-test.p-value`)),]
Huh7_cellular<-Huh7_cellular[order(-(Huh7_cellular$`DV2_Con.-Log10.Student's.T-test.p-value`)),]
Huh7_secretome<-Huh7_secretome[order(-(Huh7_secretome$`DV2_Con.-Log10.Student's.T-test.p-value`)),]

options(mc.cores=2)


network_Huh7_cellular = setRankAnalysis(Huh7_cellular$ENTREZ, collection_cellular_proteome, use.ranks=TRUE, setPCutoff = 0.01,
                                        fdrCutoff = 0.05)
exportSingleResult(network_Huh7_cellular, Huh7_cellular$ENTREZ, collection_cellular_proteome, networkName = "HuhProteome_network_ranked", IDConverter=NULL)

network_Huh7_secretome = setRankAnalysis(Huh7_secretome_foreground$ENTREZ, collection_secretome, use.ranks=TRUE, setPCutoff = 0.01,
                                         fdrCutoff = 0.05)
exportSingleResult(network_Huh7_secretome, Huh7_secretome_foreground$ENTREZ, collection_secretome, networkName = "HuhSecetome_network_ranked", IDConverter=NULL)

##Got them, but should run for all GO terms....






library(PADOG)
library(CAMERA)
library()





library(gprofiler2)
gostreCellular<-gost(
  query=Huh7_cellular_foreground$ENTREZ,
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
  custom_bg = Huh7_cellular$ENTREZ,
  numeric_ns = "ENTREZGENE_ACC",
  sources = NULL,
  as_short_link = FALSE
)

gostreSecretome<-gost(
  query=Huh7_secretome_foreground$ENTREZ,
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
  domain_scope = c("custom"),
  custom_bg = Huh7_secretome$ENTREZ,
  numeric_ns = "ENTREZGENE_ACC",
  sources = NULL,
  as_short_link = FALSE
)


Panther_run_proteome<- data.frame(Huh7_cellular$ENTREZ, Huh7_cellular$`DV2_Con.-Log10.Student's.T-test.p-value`)
Panther_run_secretome<- data.frame(Huh7_secretome$ENTREZ, Huh7_secretome$`DV2_Con.-Log10.Student's.T-test.p-value`)

write.table(Panther_run_proteome, "panther_input_proteome.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote=FALSE)
write.table(Panther_run_secretome, "panther_input_secretome.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote=FALSE)


attach(Huh7_cellular)
Huh7_cellular$fcSign=sign(`log2_DV2/Con`)

Huh7_cellular$metric=Huh7_cellular$`DV2_Con.-Log10.Student's.T-test.p-value` / Huh7_cellular$fcSign
Huh7_cellular<-Huh7_cellular[order(Huh7_cellular$metric),]
write.table(Panther_run_proteome, "panther_input_proteome.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote=FALSE)



####parameter construction
denominator<-mean(Huh7_cellular$`iDV2_Con.-Log10.Student's.T-test.p-value`)

Huh7_cellular$test_value<-(1 - Huh7_cellular$`iDV2_Con.-Log10.Student's.T-test.p-value`)/denominator




strsplit(Sys.time(), 9
