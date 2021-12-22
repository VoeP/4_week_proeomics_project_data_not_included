####Assumptions: p-values are not log10 transformed. ID:s don''t need a lot of processing, organism is human, 
####collection is generated ahead of time



library(dplyr)
library(tidyr)
library(annotate)
library(RDAVIDWebService)
library(EnrichmentBrowser)
library(gprofiler2)
library(openxlsx)
library(SetRankCorr)
library(reactome.db)
library(GeneSets.Homo.sapiens)
library(pRoloc)
uniprot2EntrezID=  createIDConverter("Homo.sapiens", "UNIPROT", "ENTREZID")

##mc.cores needs to be specified, otherwise the code won't work
options(mc.cores=2)





####Necessary functions



##accession extraction and ENTREZ conversion via wrapper

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









####Parameters
#name of table
data_name="PFP_TS-C.xlsx"
#overrepresentation or enrichment
is_ranked=FALSE
#column of gene ids
id_column="Leading_Protein"
#ranking metric for enrichment or p-value column name for overrepresentation
ranking_metric="iiTop3.tt_adjpVal.352_PFP_TS.-.352_PFP_C"
#p-value cutoff for overrepresentation
cutoff=0.05
##if results should be limited (for future)
#limit_results=c("GOMF", "GOCC", "GOBP", "REAC", "KEGG")
#limit resutls
max_nr=20

load("setr_collection_test_run.RData")

##it would be best for the user to specify the collection, since it is hard for me to test how the code reacts to
##true/false statements for the presence of the collection. Initially I ncluded the generation of the collection
##in the function, but when and if it did crash, it would always waste a lot of my time. This is especially true when
##using all the GO categories at the same time.Thus, it is advisable to either define the collection below or load a collection object.
setr_collection=
load("")
  
input_table<-read.xlsx(data_name)




####get ENTREZ

input_table=input_table %>% 
  mutate(id_column= strsplit(as.character(id_column), "-")) %>% 
  mutate(id_column= strsplit(as.character(id_column), ",")) %>% 
  unnest(id_column)
input_table$ENTREZ<-sapply(input_table[[id_column]],ProtENTREZ_bygroup,uniprot2Entrez=uniprot2EntrezID)
input_table<-input_table[which(input_table$ENTREZ != ""),]
input_table<-na.exclude(input_table)

##entirely optional: write.xlsx(input_table, file = "TS_ENTREZ_modified.xlsx")




####unranked

if (is_ranked==FALSE){
  david <- DAVIDWebService$new(email='volter.paukku@students.unibe.ch', url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  ##getting foreground and background
  foreground<-input_table[which(input_table[[ranking_metric]]<=cutoff),]
  foreground<-foreground$ENTREZ
  background<-input_table$ENTREZ
  
  ##DAVID (the email and url need to be defined each time code is run, hence they are above.)
  FG <- addList(david, foreground, idType="ENTREZ_GENE_ID", listName="isClass", listType="Gene")
  BG <- addList(david, background, idType="ENTREZ_GENE_ID", listName="all", listType="Background")
  setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
  FuncAnnotChart <- getFunctionalAnnotationChart(david)
  david_ids<-gsub("\\~.*", "", FuncAnnotChart$Term)
  david_ids<-data.frame(ids=I(david_ids), p_val=I(FuncAnnotChart$PValue))
  david_ids<- david_ids[order(david_ids$p_val),]
  david_ids<-slice_head(david_ids, n = max_nr)

  
  ##SetRank
  network<-setRankAnalysis(foreground, setr_collection, use.ranks=FALSE, setPCutoff = 0.01,
                           fdrCutoff = 0.05)
  exportSingleResult(network, foreground, setr_collection, networkName = paste(data_name, "unranked", sep = ""), IDConverter=NULL)
  SetRank_results<-read.table(paste(data_name, "unranked_pathways.txt", sep = ""), sep="\t")
  SetRank_results<-SetRank_results[-1,]
  
  ##gProfiler/GOST
  gostreUnranked<-gost(
    query=foreground,
    organism = "hsapiens",
    ordered_query = FALSE,
    multi_query = FALSE,
    significant = TRUE,
    exclude_iea = FALSE,
    measure_underrepresentation = FALSE,
    evcodes = FALSE,
    user_threshold = 0.05,
    correction_method = c("fdr"),
    domain_scope = c("annotated", "known", "custom", "custom_annotated"),
    custom_bg = input_table$ENTREZ,
    numeric_ns = "ENTREZGENE_ACC",
    sources = "GO",
    as_short_link = FALSE
  )
  gost_table<-data.frame(I(gostreUnranked[["result"]][["term_id"]]), I(gostreUnranked[["result"]][["source"]]), I(gostreUnranked[["result"]][["p_value"]]))
  gost_table<-gost_table[order(gost_table[3]),]
  gost_table<-slice_head(gost_table, n =max_nr)
  
  ##All unique ids are found and the generation of the result table is stated
  all_ids<- append((gost_table$gostreUnranked...result......term_id...), (david_ids$ids))
  all_ids<- append(all_ids, SetRank_results$V1)
  all_unique_ids<-unique(all_ids)
  result_table<-data.frame(matrix(ncol = (length(all_unique_ids))))
  ##Conversion of GO identifiers into term names
  column_names<-goIdToTerm(all_unique_ids, names=TRUE)
  colnames(result_table)<-column_names
  
  gost_pval=list()
  david_pval=list()
  setr_pval=list()
  
  ##below, the p-values for the specific IDs are found for each id in the list

  for (i in 1:length(all_unique_ids)){
    david_presence<-list()
    gost_presence<-list()
    setr_presence<-list()
    if(all_unique_ids[i] %in% david_ids$ids){
      david_presence<-filter(david_ids, ids==all_unique_ids[i])
      david_pval[i]<-david_presence$p_val

    }else{
      david_pval[i]<-"absent"
    }
    
    if(all_unique_ids[i] %in% gost_table$gostreUnranked...result......term_id...){
      gost_presence<-filter(gost_table, gostreUnranked...result......term_id...==all_unique_ids[i])
      gost_pval[i]<-gost_presence$gostreUnranked...result......p_value...

    }else{
      gost_pval[i]<-"absent"
    }
    
    if(all_unique_ids[i] %in% SetRank_results$V1){
      setr_presence<-filter(SetRank_results, V1==all_unique_ids[i])
      setr_pval[i]<-setr_presence$V7

    }else{
      setr_pval[i]<-"absent"
    }
    print("else")
  }
  
  ##Result table is completed
  result_table<-rbind(result_table, david_pval)
  result_table<-rbind(result_table, gost_pval)
  result_table<-rbind(result_table, setr_pval)
  result_table<-result_table[-1,]
  rownames(result_table)<-c("DAVID", "g:Profiler", "SetRank")
  write.table(result_table, file = "result_unranked.txt", sep = "\t")
}  














###############ranked analysis##############################

if (is_ranked==TRUE){
  input_table<-input_table[order(input_table[[ranking_metric]]),]

  ##SetRank
  network<-setRankAnalysis(input_table$ENTREZ, setr_collection, use.ranks=TRUE, setPCutoff = 0.01,
                           fdrCutoff = 0.05)
  foreground<-input_table[which(input_table[[ranking_metric]]<=cutoff),]
  exportSingleResult(network, input_table$ENTREZ, setr_collection, networkName = paste(data_name, "ranked", sep = ""), IDConverter=NULL)
  SetRank_results<-read.table(paste(data_name, "ranked_pathways.txt", sep = ""), sep="\t")
  SetRank_results<-SetRank_results[-1,]

  
  
  ##fgsea (for future)
  
  #y<-input_table[,c("ENTREZ", "metric")]
  #write.table(y,file="expression.rnk",quote=F,sep="\t",row.names=F)
  #ranks <- read.table("expression.rnk",
  #                    header=TRUE, colClasses = c("character", "numeric"))
  #ranks <- setNames(ranks$metric, ranks$ENTREZ)
  #reactome_pathways <- reactomePathways(names(ranks))
  
  #fgsea_result <- fgsea(pathways = reactome_pathways, 
  #                   stats = ranks,
  #                  minSize=15,
  #                   maxSize=500
  #)
  #fgsea_result <-fgsea_result[order(pval),]
  
  ## most significant pathways
  #significant_pathways<-fgsea_result[fgsea_result$padj<=0.01,]
  
  
  
  
  #gProfiler/GOST
  gostreRanked<-gost(
    query=input_table$ENTREZ,
    organism = "hsapiens",
    ordered_query = TRUE,
    multi_query = FALSE,
    significant = TRUE,
    exclude_iea = FALSE,
    measure_underrepresentation = FALSE,
    evcodes = FALSE,
    user_threshold = 0.05,
    correction_method = c("fdr"),
    domain_scope = c("annotated", "known", "custom", "custom_annotated"),
    custom_bg = NULL,
    numeric_ns = "ENTREZGENE_ACC",
    sources = "GO",
    as_short_link = FALSE
  )
  gost_table<-data.frame(I(gostreRanked[["result"]][["term_id"]]), I(gostreRanked[["result"]][["source"]]), I(gostreRanked[["result"]][["p_value"]]))
  gost_table<-gost_table[order(gost_table[3]),]
  gost_table<-slice_head(gost_table, n = max_nr)
  
  #All unique IDs are found, and table construction starts
  all_ids<- append(gost_table$gostreRanked...result......term_id..., SetRank_results$V1)
  all_unique_ids<-unique(all_ids)
  result_table2<-data.frame(matrix(ncol = length(all_unique_ids)))
  colnames(result_table2)<-all_unique_ids
  result_table2<-data.frame(matrix(ncol = (length(all_unique_ids)), nrow = 1))
  ##conversion of GO identifiers into term names
  column_names<-goIdToTerm(all_unique_ids, names=TRUE)
  colnames(result_table2)<-column_names

  
  ##P-values are found for each tool and ID
  gost_pval=list()
  david_pval=list()
  setr_pval=list()
  print(4)
  for (i in 1:length(all_unique_ids)){
    gost_presence<-list()
    setr_presence<-list()
    if(all_unique_ids[i] %in% gost_table$gostreRanked...result......term_id...){
      gost_presence<-filter(gost_table, gostreRanked...result......term_id...==all_unique_ids[i])
      gost_pval[i]<-gost_presence$gostreRanked...result......p_value...

    }else{
      gost_pval[i]<-"absent"
    }
    if(all_unique_ids[i] %in% SetRank_results$V1){
      setr_presence<-filter(SetRank_results, V1==all_unique_ids[i])
      setr_pval[i]<-setr_presence$V7

    }else{
      setr_pval[i]<-"absent"
    }
    

  }

  ##Final table is constructed
  result_table2<-rbind(result_table2, gost_pval)
  result_table2<-rbind(result_table2, setr_pval)
  result_table2<-result_table2[-1,]
  rownames(result_table2)<-c("g:Profiler", "SetRank")
  write.table(result_table2, file = "result_ranked.txt", sep = "\t")
  
}


