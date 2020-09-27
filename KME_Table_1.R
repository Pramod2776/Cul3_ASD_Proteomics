library(tidyverse)
library(limma)
library(edgeR)
library(UniProt.ws)
library(WGCNA)
library(rentrez)
library("biomaRt")
library(doParallel)
enableWGCNAThreads()
allowWGCNAThreads()
registerDoParallel()
options(tibble.width=Inf)

RData_files <- list.files(path = "./Input_files", full.names = T)
ms <- UniProt.ws::UniProt.ws(taxId = 10090)

for (fn in RData_files){
  comp_file_name <- str_split(str_split(fn, "/")[[1]][3], "_")
  Period <- comp_file_name[[1]][1]
  Region <- comp_file_name[[1]][2]
  
  load(fn)
  
  geneModuleMembership<- signedKME(t(networks$datExpr), datME = networks$MEs$eigengenes,corFnc = "bicor")
  labels2colors(as.numeric(gsub("ME", "", colnames(networks$MEs$eigengenes))))
  colnames(geneModuleMembership)=paste("kME","_", labels2colors(as.numeric(gsub("ME", "", colnames(networks$MEs$eigengenes)))),sep="")
  genes<- rownames(networks$datExpr)
  modules <- networks$merged$colors
  module_number <- as.numeric(gsub("ME", "", colnames(networks$MEs$eigengenes)))
  module_number <- module_number %>%
    set_names(as.character(module_number))
  gene_module<- lapply(set_names(names(module_number)),function(n){
    data.frame(
      Gene_name = genes[modules==module_number[[n]]],
      Module = paste0(as.numeric(n),labels2colors(module_number[[n]]))
    )
  }
  )%>%
    bind_rows()
  
  proteins <- gene_module$Gene_name
  proteins1<-split(proteins, ceiling(seq_along(proteins)/100))
  columns <- c("ENTRY-NAME", "PROTEIN-NAMES", "GENES", "ENSEMBL", "ENSEMBL_PROTEIN")
  
  mouse_anno2<-(lapply(names(proteins1), function(names){
    mouse_anno <- UniProt.ws::select(ms, proteins1[[names]], columns, "UNIPROTKB") %>%
      mutate(GeneSymbol = gsub(" .*", "", GENES)) %>%
      dplyr::select(
        Accessions = "UNIPROTKB", Gene = "ENTRY-NAME",
        Description = "PROTEIN-NAMES", ENSEMBL= "ENSEMBL", ENSEMBL_PROTEIN="ENSEMBL_PROTEIN", GeneSymbol
      )
    
  })) %>%
    bind_rows()
  
  gene_module_ensembl<- mouse_anno2 %>%
    distinct(Gene, .keep_all = TRUE)
  
  gene_module <- gene_module %>%
    rename(replace = c("Gene_name" = "Accessions"))
  gene_module_ensembl_Acc <- gene_module %>%
    left_join(gene_module_ensembl, by = "Accessions")%>%
    dplyr::select(c("Accessions", "GeneSymbol", "ENSEMBL", "Module"))
  
  sort_genemodule_membership<- geneModuleMembership[match(gene_module_ensembl_Acc$Accessions,(rownames(geneModuleMembership))),]
  genemodule_membership_KME<- cbind(gene_module_ensembl_Acc, sort_genemodule_membership)
  
  genemodule_membership_KME_table <- genemodule_membership_KME%>%
    arrange(GeneSymbol)%>%
    arrange(Module)
  
  writexl::write_xlsx(genemodule_membership_KME_table, 
                      paste0(Period,"_",Region,"_","WGCNA_Proteomics_KME_Table.xlsx"),
                      col_names = T)
  
}

pre_KME_fns <- list.files(".","WGCNA")
ensembl <- useEnsembl(biomart = "ensembl", dataset = 'mmusculus_gene_ensembl')

for (fn in pre_KME_fns){
  input_genes <- readxl::read_xlsx(fn, sheet = 1) %>% 
    pull(GeneSymbol)
  
  Ensembl_IDs <- getBM(attributes=c('mgi_symbol', 'ensembl_gene_id'), 
        filters = 'mgi_symbol', 
        values = input_genes, 
        mart = ensembl) %>%
    rename(replace = c('mgi_symbol' = "GeneSymbol"))
  
  merged_sheet <- left_join(readxl::read_xlsx(fn, sheet = 1), 
                                    Ensembl_IDs, 
                                    by = "GeneSymbol")
  
  merged_sheet$ENSEMBL[is.na(merged_sheet$ENSEMBL)] <-  
    merged_sheet$ensembl_gene_id[is.na(merged_sheet$ENSEMBL)]
  
  new_KME_table <- merged_sheet[!duplicated(paste(merged_sheet$Accessions,
                                                  merged_sheet$GeneSymbol,
                                                  merged_sheet$Module)),]
  new_KME_table$ensembl_gene_id = NULL
  
  remove_na_table <- new_KME_table[!is.na(new_KME_table$ENSEMBL),]
  
  unfilled_ensembl_list <- new_KME_table[is.na(new_KME_table$ENSEMBL),]
  
  writexl::write_xlsx(new_KME_table, fn, col_names = T)
  writexl::write_xlsx(remove_na_table, paste0("./Output_files/NA_removed_",fn), col_names = T)
  writexl::write_xlsx(unfilled_ensembl_list, paste0("./Output_files/Unfilled_info_",fn), col_names = T)
  
  
}











