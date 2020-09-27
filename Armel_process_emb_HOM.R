setwd("/Users/pramod/Desktop/Proteomics/Fwd__SL_combat_script_validate/Cul3 - regulated PG tables/Embryonic//")
library(tidyverse)
library(doParallel)
registerDoParallel()
EMHOM<-readxl::read_xlsx("reg. proteinGroups - JM-Embryonic-2 HOM.xlsx", sheet = 1)
EMHOM_list<-list()
counter<-1
labels <- c(10, 20, 30)
for (sl in c(293,779,944)){
  KO <- EMHOM %>%
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
  EMHOM_list[[counter]] <- KO
  EMHOM_list[[counter + 1]] <- unique_KO
  counter <- counter + 2
  
}

names(EMHOM_list) <- sapply(
  labels, 
  function(l) paste(c("EMHOM_list", "EMHOM_unique_data"), l, sep = "_")
) %>%
  as.character()
names(EMHOM_list$EMHOM_unique_data_10)<-c("Genes_10")
names(EMHOM_list$EMHOM_unique_data_20)<-c("Genes_20")
names(EMHOM_list$EMHOM_unique_data_30)<-c("Genes_30")


writexl::write_xlsx(EMHOM_list$EMHOM_unique_data_10, "reg_proteinGroups-JM-Embryonic-2_HOM_10_mod.xlsx", col_names = T)
writexl::write_xlsx(EMHOM_list$EMHOM_unique_data_20, "reg_proteinGroups-JM-Embryonic-2_HOM_20_mod.xlsx", col_names = T)
writexl::write_xlsx(EMHOM_list$EMHOM_unique_data_30, "reg_proteinGroups-JM-Embryonic-2_HOM_30_mod.xlsx", col_names = T)

writexl::write_xlsx(EMHOM_list$EMHOM_list_10, "reg_proteinGroups-JM-Embryonic-2_HOM_unique_list_10_mod_A.xlsx", col_names = T)
writexl::write_xlsx(EMHOM_list$EMHOM_list_20, "reg_proteinGroups-JM-Embryonic-2_HOM_unique_list_20_mod_A.xlsx", col_names = T)



