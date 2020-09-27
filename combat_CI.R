setwd("/Users/pramod//Desktop/Proteomics/Fwd__SL_combat_script_validate/CI_Intensity/")
install.packages("BiocManager")
BiocManager::install("tidyverse")
library(tidyverse)
BiocManager::install("writexl")
library(writexl)
BiocManager::install("qPLEXanalyzer")
library(qPLEXanalyzer)
BiocManager::install("gridExtra")
library(gridExtra)
BiocManager::install("UniProt.ws")
library(UniProt.ws)
BiocManager::install("dplyr")
library(dplyr)
BiocManager::install("edgeR")
library(edgeR) 
BiocManager::install("sva")
library(sva)
BiocManager::install("limma")
library(limma) 


dat<-read.table("census-out-Cul3_E17.5-16719_cortex.txt", header = T,  sep = "\t", 
                fill = TRUE)

filter_protein <- dat %>%
  dplyr::select(PROTEIN, SEQUENCE) %>%
  distinct(PROTEIN, SEQUENCE) %>%
  count(PROTEIN) %>%
  filter(n >= 2) %>%
  pull(PROTEIN) %>%
  as.character()

dat_filt <- dat%>%
  filter(PROTEIN %in% filter_protein) %>%
  dplyr::select(-starts_with("m"))%>%
  dplyr::select(
    PROTEIN, SEQUENCE, UNIQUE,starts_with("norm")
  ) %>%
  mutate(PROTEIN = as.character(PROTEIN))
names(dat_filt)<-c("PROTEIN", "SEQUENCE", "UNIQUE", "126_S2", "127N_S2", "127C_S4", 
              "128N_S4", "128C_S4", "129N_S4", "129C_S3", "130N_S3", "130C", "131N")


dat_filt<-dat_filt%>%
  filter(! str_detect(PROTEIN, "contaminant|Reverse"))

dim(dat_filt)
meta<-read.csv("cul3_cortex_meta.csv", header = T)
meta<-meta[order(meta$SampleGroup),]
meta<-meta %>%
  mutate(Sample=c("127N_S2", "128N_S4", "129N_S4", "130N_S3", "130C", "131N", "126_S2","127C_S4", 
                  "128C_S4", "129C_S3"))
dat_filt_match<-dat_filt[,4:13]
dat_filt_match<-dat_filt_match[match(meta$Sample, colnames(dat_filt_match))]
dat_sort<-cbind(dat_filt[,1:3], dat_filt_match)
dat_sort_A<-dat_sort[, c("PROTEIN", "SEQUENCE", "UNIQUE", "127N_S2", "128N_S4", "129N_S4", "130N_S3", "126_S2","127C_S4", 
                         "128C_S4", "129C_S3","130C", "131N") ]
meta_group<-read.csv("meta_group.csv", header = T)
meta_group<-meta_group %>%
  mutate(Sample=c("127N_S2", "128N_S4", "129N_S4", "130N_S3", "126_S2","127C_S4", 
                  "128C_S4", "129C_S3", "130C", "131N"))
meta_group<-meta_group %>%
  mutate(SampleName=Sample)
exp2=list(intensities=dat_sort_A, metadata=meta_group)
Mnset<-convertToMSnset(exp2$intensities, metadata = exp2$metadata, indExpData = c(4:13), 
                       Sequences = 2, Accessions = 1)
exprs(Mnset)<-exprs(Mnset)+0.01

proteins <- unique(fData(Mnset)$Accessions)
columns <- c("ENTRY-NAME", "PROTEIN-NAMES", "GENES")
ms <- UniProt.ws::UniProt.ws(taxId = 10090)
mouse_anno <- UniProt.ws::select(ms, proteins, columns, "UNIPROTKB") %>%
  mutate(GeneSymbol = gsub(" .*", "", GENES)) %>%
  dplyr::select(
    Accessions = "UNIPROTKB", Gene = "ENTRY-NAME",
    Description = "PROTEIN-NAMES", GeneSymbol
  )

MSnset_Sum <- summarizeIntensities(Mnset, sum, mouse_anno)
unorm_sum_intensity<-exprs(MSnset_Sum)

data_start <- unorm_sum_intensity
data_start<-data_start[, 1:8]

data_raw <- na.omit(data_start)
colnames<-c( "Cul3_HET_1", "Cul3_HET_2", "Cul3_HET_3", "Cul3_HET_4", "WT_1", "WT_2", 
             "WT_3", "WT_4")
names(data_raw)<-colnames
data_raw<- as.data.frame(data_raw)
meta_PCA<- meta_group[-c(9:10),]
prc <- prcomp(t(data_raw), scale = TRUE, center = TRUE)
prc$x %>%
  as.data.frame() %>%
  mutate(SampleName = rownames(.)) %>%
  inner_join(meta_PCA, by = "SampleName") %>%
  separate(SampleID, into = c("S", "SID"), sep = "-") %>%
  ggplot(
    data = .,
    aes(
      x = PC1, y = PC2, colour = SampleGroup, shape = S, label = SampleName
    )
  ) +
  geom_text(vjust = 1) +
  geom_point()

cul3_raw<-data_raw["Q9JLV5",] %>%
  reshape2::melt() %>%
  mutate(
    Genotype = ifelse(
      str_detect(variable, "N"),
      "HET", "WT"
    )
  ) %>%
  ggplot(
    data = .,
    aes(
      x = variable, y = value, fill = Genotype
    )
  ) +
  facet_grid(. ~ Genotype, scales = "free_x") +
  geom_bar(stat = "identity", colour = "black") +
  ggtitle("normalized intensity")

cul3_raw

design <- meta_group[-c(9:10),]
design$SampleGroup <- factor(design$SampleGroup, levels = c("Cul3_HET", "WT"))
mod <- model.matrix(~ SampleGroup + Gender, data = design)
batch <- design$Batch
data_raw=as.matrix(data_raw)
data_combat <- ComBat(dat = data_raw, batch = batch, mod = mod, par.prior = TRUE)
data_combat <-  data_combat[apply(data_combat, 1, function(x) all(x > 0)), ]
data_combat_A<-as.data.frame(data_combat)
cul3_combat<-data_combat_A["Q9JLV5",] %>%
  reshape2::melt() %>%
  mutate(
    Genotype = ifelse(
      str_detect(variable, "N"),
      "HET", "WT"
    )
  ) %>%
  ggplot(
    data = .,
    aes(
      x = variable, y = value, fill = Genotype
    )
  ) +
  facet_grid(. ~ Genotype, scales = "free_x") +
  geom_bar(stat = "identity", colour = "black") +
  ggtitle("Combat normalized on norm intensities")

cul3_combat
boxplot(log2(data_combat), notch = TRUE, col = rep(c("red", "green"), each = 4), 
        main = "ComBat batch correction of CI data\nExp1 (red), Exp2 (green)")

plotDensities(log2(data_combat), col = rep(c("red", "green"), each=4), main = "ComBat CI data")
plotMDS(log2(data_combat), main = "normalized Combat CI data with Gender")

data_combat_pca<-data_combat
prc <- prcomp(t(data_combat_pca), scale = TRUE, center = TRUE)
prc$x %>%
  as.data.frame() %>%
  mutate(SampleName = rownames(.)) %>%
  inner_join(meta_PCA, by = "SampleName") %>%
  separate(SampleID, into = c("S", "SID"), sep = "-") %>%
  ggplot(
    data = .,
    aes(
      x = PC1, y = PC2, colour = SampleGroup, shape = S, label = SampleName
    )
  ) +
  geom_text(vjust = 1) +
  geom_point()

group <- rep(c("CUL3_HET", "WT"), each=4)
group <- factor(group, levels = c("WT", "CUL3_HET"))
y_sl <- DGEList(counts = data_combat, group = group)
y_sl$samples
y_sl <- estimateDisp(y_sl)
plotBCV(y_sl)
et_sl <- exactTest(y_sl, pair = c("WT", "CUL3_HET"))
summary(decideTestsDGE(et_sl))
tt_sl <- topTags(et_sl, n = Inf, sort.by = "none")
tt_sl <- tt_sl$table
dim(tt_sl)
ggplot(tt_sl, aes(PValue))+
  geom_histogram(bins=100, fill="white", color="black")+
  geom_hline(yintercept = mean(hist(tt_sl$PValue, breaks = 100, plot = FALSE)$counts[26:100]))+
  ggtitle("p-value distribution CI intensity edgeR with Gender")
mouse_anno_sort=mouse_anno[match(rownames(tt_sl), mouse_anno$Accessions),]
tt_sl_anno=cbind(tt_sl, mouse_anno_sort)
tt_sl_anno=tt_sl_anno%>%
  dplyr::arrange(FDR)%>%
  dplyr::select(Accessions, Gene, GeneSymbol, Description, logFC, logCPM, PValue, FDR)%>%
  arrange(FDR)
DPE=tt_sl_anno[tt_sl_anno$FDR < 0.1,]
nrow(DPE)
tail(DPE)
write.csv(tt_sl_anno, "CI_intensity_DPE_cul3_anno_order_gender_combat_A.csv",)

mod_mat <- mod
colnames(mod_mat) <- make.names(colnames(mod_mat))

lm_fit <- lmFit(log2(data_combat), design = mod_mat)
contr_mat <- makeContrasts(
  HET = -SampleGroupWT,
  levels = colnames(mod_mat)
)
contr_fit <- contrasts.fit(lm_fit, contr_mat)
ebayes <- eBayes(contr_fit)
tt <- topTable(ebayes, coef = "HET", number = Inf)
dpe <- tt %>%
  filter(adj.P.Val <= 0.05 & abs(logFC) >= log2(1.5))
ggplot(tt, aes(P.Value))+
  geom_histogram(bins=100, fill="white", color="black")+
  geom_hline(yintercept = mean(hist(tt$P.Value, breaks = 100, plot = FALSE)$counts[26:100]))+
  ggtitle("p-value distribution CI intensity limma with gender")
mouse_anno_sort_limma=mouse_anno[match(rownames(tt), mouse_anno$Accessions),]
tt_limma<-cbind(tt, mouse_anno_sort_limma)
DPE_limma=tt_limma[tt_limma$adj.P.Val < 0.1,]
nrow(DPE_limma)
head(DPE_limma)
write.csv(tt_limma, "CI_intensity_limma_result_filter_make_contrasts_gender_combat_anno_A.csv",)


volcano_input <- bind_rows(
  tt_limma %>%
    mutate(Protein = rownames(.)) %>%
    dplyr::rename(fdr = adj.P.Val) %>%
    mutate(Method = "LIMMA") %>%
    dplyr::select(Protein, logFC, fdr, Method) %>% 
    mutate(Upregulated = (tt_limma$logFC > 0.58 & tt_limma$adj.P.Val<0.05)) %>%
    mutate(Downregulated = (tt_limma$logFC < -0.58 & tt_limma$adj.P.Val<0.05)) %>%
    mutate(Gene=tt_limma$GeneSymbol),
  tt_sl_anno %>%
    mutate(Protein = rownames(.)) %>%
    dplyr::rename(fdr = FDR) %>%
    mutate(Method = "EDGER") %>%
    dplyr::select(Protein, logFC, fdr, Method) %>% 
    mutate(Upregulated = (tt_sl_anno$logFC > 0.58 & tt_sl_anno$FDR<0.05)) %>%
    mutate(Downregulated = (tt_sl_anno$logFC < -0.58 & tt_sl_anno$FDR<0.05)) %>%
    mutate(Gene=tt_sl_anno$GeneSymbol)
) %>%
  mutate(Direction = mapply(
    function(up, down) {
      if (up) {
        return("UP")
      } else if (down) {
        return("DOWN")
      } else {
        return("NA")
      }
    },
    Upregulated,
    Downregulated
  )) %>%
  filter(! str_detect(Protein, "contaminant|Reverse"))

ggplot() +
  facet_wrap(~ Method) +
  geom_point(
    data = volcano_input %>%
      filter(!Upregulated & !Downregulated),
    aes(
      x = logFC, y = -log10(fdr)
    )
  ) +
  geom_point(
    data = volcano_input %>%
      filter(Upregulated | Downregulated),
    aes(
      x = logFC, y = -log10(fdr), colour = Direction
    )
  ) +
  scale_colour_manual(
    values = c(
      "UP" = "red",
      "DOWN" = "blue",
      "NA" = "lightgrey"
    )
  ) +
  ggtitle(label = "Cul3vsWT") +
  theme_bw(base_size = 16)+
  geom_vline(xintercept = 0.58, colour = "black")+
  geom_vline(xintercept = -0.58, colour = "black")+
  geom_hline(yintercept = 1.3, colour = "red")+
  scale_y_continuous(trans = "log1p")+
  ggrepel::geom_text_repel(
    data = volcano_input %>%
      filter(Upregulated | Downregulated), 
    aes(label = Gene, x = logFC, y = -log10(fdr))
  )

dpe_full_merged <- merge(
  tt %>%
    mutate(Protein = rownames(.)),
  tt_sl %>%
    mutate(Protein = rownames(.)),
  by = "Protein",
  suffixes = c(".LIMMA", ".EDGER")
)

ggplot(
  data = dpe_full_merged,
  aes(
    x = logFC.LIMMA, y = logFC.EDGER
  )
) +
  geom_point()

correlation<-cor(dpe_full_merged$logFC.LIMMA, dpe_full_merged$logFC.EDGER, method = "kendall")
correlation_test<-cor.test(dpe_full_merged$logFC.LIMMA, dpe_full_merged$logFC.EDGER, method = "kendall")
















