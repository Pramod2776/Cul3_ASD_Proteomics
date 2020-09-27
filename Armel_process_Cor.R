setwd("/Users/pramod/Desktop/Proteomics/Fwd__SL_combat_script_validate/Cul3 - regulated PG tables/Adult/")
library(tidyverse)
library(doParallel)
registerDoParallel()
ACor<-readxl::read_xlsx("reg. proteinGroups - KO Cortex.xlsx", sheet = 1)
ACor_data <- list()
counter <- 1
labels <- c("10", "20", "30")
for (sl in c(96, 133, 135)) {
  Cor <- ACor %>%
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
  unique_Cor <- Cor %>%
    filter(!is.na(Genes)) %>%
    pull(Genes) %>%
    unique()%>%
    as.data.frame()
  ACor_data[[counter]] <- Cor
  ACor_data[[counter + 1]] <- unique_Cor
  counter <- counter + 2
}
names(ACor_data) <- sapply(
  labels, 
  function(l) paste(c("ACor_data", "ACor_unique_data"), l, sep = "_")
) %>%
  as.character()

ACor_data

names(ACor_data$ACor_unique_data_10)<-c("Genes_10")
names(ACor_data$ACor_unique_data_20)<-c("Genes_20")
names(ACor_data$ACor_unique_data_30)<-c("Genes_30")



writexl::write_xlsx(ACor_data$ACor_unique_data_10, "reg_proteinGroups-KO_Adult_sorted_Cortex_10_mod.xlsx", col_names = T)
writexl::write_xlsx(ACor_data$ACor_unique_data_20, "reg_proteinGroups-KO_Adult_sorted_Cortex_20_mod.xlsx", col_names = T)
writexl::write_xlsx(ACor_data$ACor_unique_data_30, "reg_proteinGroups-KO_Adult_sorted_Cortex_30_mod.xlsx", col_names = T)

writexl::write_xlsx(ACor_data, "reg_proteinGroups-KO_Adult_sorted_Cortex_mod_A.xlsx", col_names = T)

 
