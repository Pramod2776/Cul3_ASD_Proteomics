setwd("/Users/pramod/Desktop/Proteomics/Fwd__SL_combat_script_validate/ASD_enrichment/")
library(tidyverse)
Compare_list<-lapply(
  setNames(
    nm = readxl::excel_sheets("ListsFromLiterature.xlsx")
  ), 
  function(sheet){
    readxl::read_xlsx("ListsFromLiterature.xlsx",
                      sheet = sheet)
  }
)


DE<-readxl::read_xlsx("CI_SL_intensity_limma_result_filter_make_contrasts_gender_combat_anno_A_filter.xlsx", 
                  sheet = 1)
DPE<-as.character(DE %>%
  filter(adj.P.Val < 0.1)%>%
  pull(GeneSymbol))

NDPE<-as.character(DE %>%
                    filter(adj.P.Val >= 0.1)%>%
                    pull(GeneSymbol))

names<-names(Compare_list)
FT<-lapply(
  names,
  function(ncl){
    cl<-na.omit(Compare_list[[ncl]][[1]])
    input=na.omit(toupper(c(DPE, NDPE)))
    bg=unique(c(input, cl))
    a=length(intersect(toupper(DPE), cl))
    c=length(intersect(DPE, bg[!bg %in% cl]))
    b=length(intersect(bg[!bg %in% toupper(DPE)], cl))
    d=length(bg[!bg %in% c(DPE, cl)])
    contingency= matrix(
      c(a,b,c,d), nrow = 2, ncol = 2, byrow = TRUE
    )
    result=fisher.test(contingency)%>%
      broom::tidy()%>%
      mutate(List=ncl)
  
  }
) %>%
  bind_rows()%>%
  mutate(FDR=p.adjust(p.value, method = "fdr"))
  
writexl::write_xlsx(FT, "Gene_list_enrichment_E17.5_CX_mod.xlsx", col_names = T)


















