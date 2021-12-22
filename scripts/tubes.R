##the remaining sets, tube and carrier

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


PFP_TS_foreground<-subset(PFP_TS, `iiTop3.tt_adjpVal.352_PFP_TS.-.352_PFP_C`<0.05)


PFP_C$ProtsACs=sapply(PFP_C[,"Protein.group"],ExtractACs)
PFP_C<-PFP_C[!is.na(PFP_C[,"log.ANOVA.p.value"]),]

PFP_C=PFP_C %>% 
  mutate(Protein.group= strsplit(as.character(Protein.group), ";")) %>% 
  unnest(Protein.group)
PFP_C$ENTREZ<-sapply(PFP_C$Protein.group,ProtENTREZ_bygroup,uniprot2Entrez=uniprot2EntrezID)
PFP_C_noENTREZ<-PFP_C[which(PFP_C$ENTREZ==""),]
PFP_C<-PFP_C[which(PFP_C$ENTREZ != ""),]
PFP_C_foreground<-subset(PFP_C, ANOVA.p.value<0.05)





PFP_C_foreground<-PFP_C_foreground[order((PFP_C_foreground$ANOVA.p.value)),]
PFP_TS_foreground<-PFP_TS_foreground[order((PFP_TS_foreground$`iiTop3.tt_adjpVal.352_PFP_TS.-.352_PFP_C`)),]


Panther_run_carrier<- data.frame(PFP_C_foreground$ENTREZ, PFP_C_foreground$ANOVA.p.value)
Panther_run_tube<- data.frame(PFP_TS_foreground$ENTREZ, PFP_TS_foreground$`iiTop3.tt_adjpVal.352_PFP_TS.-.352_PFP_C`)

write.table(Panther_run_carrier, "panther_input_carrier.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote=FALSE)
write.table(Panther_run_tube, "panther_input_tube.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote=FALSE)


load("collection_carrier_500_GOBP.RData")
load("collection_tube_500_GOBP.RData")


options(mc.cores=2)
library(SetRankCorr)

network_ranked_tube = setRankAnalysis(unlistPFP_TS_foreground$ENTREZ, collection_tube, use.ranks=TRUE, setPCutoff = 0.01,
                                        fdrCutoff = 0.05)

exportSingleResult(network_ranked_tube, PFP_TS_foreground, collection_tube, networkName = "tube_system_network_ranked", IDConverter=NULL)



network_ranked_tube = setRankAnalysis(PFP_C_foreground$ENTREZ, collection_carrier, use.ranks=TRUE, setPCutoff = 0.01,
                                        fdrCutoff = 0.05)

exportSingleResult(network_ranked_tube, PFP_TS_foreground, collection_tube, networkName = "tube_system_network_ranked", IDConverter=NULL)

library(CAMERA)
