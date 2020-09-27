setwd("/Users/pramod/Desktop/Proteomics/Fwd__SL_combat_script_validate/simon/")
library(tidyverse)
library(ggrepel)
library(ggplot2)
data<-readxl::read_xlsx("CI_SL_intensity_limma_result_filter_make_contrasts_gender_combat_anno_A_filter.xlsx", sheet = 1)
Vol_data<-data %>%
  dplyr::select(Protein_ID,  adj.P.Val, logFC, GeneSymbol)%>%
  dplyr::rename(Protein_ID = Protein_ID, pvalue = adj.P.Val, lfc = logFC, gene = GeneSymbol)
GO_term<-readxl::read_xlsx("C6-E17.5 CX-Proteomics-15-ALL_Poster-2.xlsx", sheet = 2) 
neuron<-GO_term %>%
  transform(intersection=strsplit(intersection, ","))%>%
  unnest()
neuron_neupro<-neuron%>%
  filter(GO_Term=="neuron projection development")%>%
  pull(intersection)

neuron_interfil<-neuron%>%
  filter(GO_Term=="intermediate filament cytoskeleton organization")%>%
  pull(intersection)

neuron_regtrans<-neuron%>%
  filter(GO_Term=="regulation of trans-synaptic signaling")%>%
  pull(intersection)

neuron_actfil<-neuron%>%
  filter(GO_Term=="actin filament polymerization")%>%
  pull(intersection)

genes16p_down<-Vol_data%>%
  filter(lfc <= 0)%>%
  filter(toupper(gene) %in% neuron_neupro)

genes16p_up<-Vol_data%>%
  filter(lfc >= 0)%>%
  filter(toupper(gene) %in% neuron_neupro)

label<-c("MEF2", "NFIB", "NEFL", "VIM", "STXBP1", "LIMK1", "STX1B", "MECP2",
         "PLXNB2", "ARHGAP33", "TRIM67", "ARF6", "MAP1B", "PLXNB1", "MAP6", 
         "TUBB3", "ANK3", "SYNE1", "CUL3")

label_neuron_neupro_down<-Vol_data%>%
  filter(lfc <= 0)%>%
  filter(toupper(gene) %in% label)%>%
  pull(gene)

label_neuron_neupro_up<-Vol_data%>%
  filter(lfc >= 0)%>%
  filter(toupper(gene) %in% label)%>%
  pull(gene)

Vol_data$pvalue=-log10(Vol_data$pvalue)

volcano_DUPvsDEL <- Vol_data %>% mutate(color = ifelse (Vol_data$gene %in% label_neuron_neupro_up,
                                                        yes = "label_up",
                                                        no = ifelse (Vol_data$gene %in% label_neuron_neupro_down,
                                                                     yes ="label_down",
                                                        no = ifelse (Vol_data$gene %in% genes16p_down$gene,
                                                                     yes ="genes16p_down",
                                                                     no = ifelse (Vol_data$gene %in% genes16p_up$gene,
                                                                                  yes ="genes16p_up",
                                                                                  no = ifelse(Vol_data$lfc > 0.38 & Vol_data$pvalue > 1.0,
                                                                                              yes = "Upregulated", 
                                                                                              no = ifelse(Vol_data$lfc < -0.38 & Vol_data$pvalue > 1.0, 
                                                                                                          yes = "Downregulated",
                                                                                                          no="none")))))))




volcano_DUPvsDEL[volcano_DUPvsDEL$color == "label_down",]

volcano_DUPvsDEL_none<-volcano_DUPvsDEL[volcano_DUPvsDEL$color == "none",]
volcano_DUPvsDEL_up<-volcano_DUPvsDEL[volcano_DUPvsDEL$color == "Upregulated",]
volcano_DUPvsDEL_down<-volcano_DUPvsDEL[volcano_DUPvsDEL$color == "Downregulated",]
volcano_DUPvsDEL_16p_down<-volcano_DUPvsDEL[volcano_DUPvsDEL$color == "genes16p_down",]
volcano_DUPvsDEL_label_down<-volcano_DUPvsDEL[volcano_DUPvsDEL$color == "label_down",]
volcano_DUPvsDEL_label_up<-volcano_DUPvsDEL[volcano_DUPvsDEL$color == "label_up",]
volcano_DUPvsDEL_16p_up<-volcano_DUPvsDEL[volcano_DUPvsDEL$color == "genes16p_up",]


pdf("Vol+plot_CX_E17.5_proteomics_simon_neuro_proj_cul3.pdf", width = 16, height = 20, useDingbats = FALSE)
  
  ggplot() + 
  geom_point(data = volcano_DUPvsDEL_none, aes(x = lfc, y = pvalue, color = factor(color)), size = 4, alpha = 0.4, na.rm = T) +
  geom_point(data = volcano_DUPvsDEL_up, aes(x = lfc, y = pvalue, color = factor(color)), size = 4, alpha = 0.4, na.rm = T)+
  geom_point(data = volcano_DUPvsDEL_down, aes(x = lfc, y = pvalue, color = factor(color)), size = 4, alpha = 0.4, na.rm = T)+
  geom_point(data = volcano_DUPvsDEL_16p_down, aes(x = lfc, y = pvalue, color = factor(color)), size = 7, alpha = 0.8, na.rm = T)+
  geom_point(data = volcano_DUPvsDEL_16p_up, aes(x = lfc, y = pvalue, color = factor(color)), size = 7, alpha = 0.8, na.rm = T)+
  geom_point(data = volcano_DUPvsDEL_label_down, aes(x = lfc, y = pvalue, color = factor(color)), size = 7, alpha = 0.8, na.rm = T)+
  geom_point(data = volcano_DUPvsDEL_label_up, aes(x = lfc, y = pvalue, color = factor(color)), size = 7, alpha = 0.8, na.rm = T)+
  
  
  # add gene points
  theme_bw(base_size = 64) + # clean up theme
  theme(legend.position = "none") + # remove legend 
  ggtitle(label = "CX_E17.5_Neuro_Proj") +  # add title
  #xlab(expression(log[2]("DUP" / "CTL"))) + # x-axis label
  xlab("log2 FC (HET/WT)")+
  ylab(expression(-log[10]("FDR"))) + # y-axis label
  geom_vline(xintercept = 0.38, colour = "black") + # add line at 0
  geom_vline(xintercept = -0.38, colour = "black") + # add line at 0
  geom_hline(yintercept = 1.3, colour = "black") + # p(0.01) =2
  scale_y_continuous(trans = "log1p") +
  scale_color_manual(values = c("Upregulated" = "grey60", 
                                "Downregulated" ="grey60",
                                "none" = "grey60", "genes16p_down"="#377EB8", "genes16p_up"="#E41A1C", "label_down"="#377EB8", "label_up"="#E41A1C"))+
  geom_text_repel(
    # data=subset(volcano_DUPvsDEL,gene_id %in% genes16psubset),
    data = subset(volcano_DUPvsDEL, color %in% c( "label_down", "label_up")), 
    aes(label = gene, x = lfc, y = pvalue),
    col="black",size=10.5,fontface="bold",box.padding = unit(1.5, "lines"),
    point.padding = unit(0.5, "lines"))+ geom_text(check_overlap = TRUE)

dev.off()


