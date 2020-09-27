setwd("/Users/pramod/Desktop/Proteomics/Fwd__SL_combat_script_validate/Cul3 - regulated PG tables/Embryonic//")
library(tidyverse)
library(doParallel)
registerDoParallel()
EMKO<-readxl::read_xlsx("reg. proteinGroups - JM-Embryonic-1 KO.xlsx", sheet = 1)
EMKO_list<-list()
counter<-1
labels <- c(10, 20, 30)
for (sl in c(166,179,184)){
  KO <- EMKO %>%
    slice(1:sl)%>%
    transform(Genes = strsplit(as.character(Genes), ";"))%>%
    unnest(Genes)%>%
    transform(RazorPeptideIDs = strsplit(as.character(RazorPeptideIDs), ";"))%>%
    unnest(RazorPeptideIDs)%>%
    transform(ProteinIDs = strsplit(as.character(ProteinIDs), ";"))%>%
    unnest(ProteinIDs)%>%
    transform(UniquePeptideIDs = strsplit(as.character(UniquePeptideIDs), ";"))%>%
    unnest(UniquePeptideIDs)%>%
    dplyr::select(Common.Names, ProteinIDs, Genes, PeptideIDs, RazorPeptideIDs, UniquePeptideIDs, everything()
    )
  unique_KO <- KO %>%
    filter(!is.na(Genes)) %>%
    pull(Genes) %>%
    unique()%>%
    as.data.frame()
  EMKO_list[[counter]] <- KO
  EMKO_list[[counter + 1]] <- unique_KO
  counter <- counter + 2
  
}

names(EMKO_list) <- sapply(
  labels, 
  function(l) paste(c("EMKO_list", "EMKO_unique_data"), l, sep = "_")
) %>%
  as.character()
names(EMKO_list$EMKO_unique_data_10)<-c("Genes_10")
names(EMKO_list$EMKO_unique_data_20)<-c("Genes_20")
names(EMKO_list$EMKO_unique_data_30)<-c("Genes_30")


writexl::write_xlsx(EMKO_list$EMKO_unique_data_10, "reg_proteinGroups-JM-Embryonic-1_KO_10_mod.xlsx", col_names = T)
writexl::write_xlsx(EMKO_list$EMKO_unique_data_20, "reg_proteinGroups-JM-Embryonic-1_KO_20_mod.xlsx", col_names = T)
writexl::write_xlsx(EMKO_list$EMKO_unique_data_30, "reg_proteinGroups-JM-Embryonic-1_KO_30_mod.xlsx", col_names = T)

writexl::write_xlsx(EMKO_list, "reg_proteinGroups-JM-Embryonic-1_KO_mod_A.xlsx", col_names = T)



