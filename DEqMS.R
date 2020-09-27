setwd("/Users/pramod//Desktop/Proteomics/Fwd__SL_combat_script_validate/SL_norm/")
library(tidyverse)
library(writexl)
library(qPLEXanalyzer)
library(gridExtra)
library(UniProt.ws)
library(dplyr)
library(edgeR) 
library(sva)
library(limma) 
library(doParallel)
registerDoParallel()
BiocManager::install("DEqMS")
library(DEqMS)

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

#psm count

psm.count.table = as.data.frame(table(dat_filt$PROTEIN))
rownames(psm.count.table)=psm.count.table$Var1
write.csv(psm.count.table, "psm.count.table.csv",)

plot(sort(log2(psm.count.table$Freq)),pch=".",
     xlab="Proteins ordered by PSM count",ylab="log2(PSM count)")

# preprocessing

dat_filt%>%
  count(PROTEIN)%>%
  filter(PROTEIN == "Q9JLV5")
dat_filt%>%
  count(PROTEIN)%>%
  filter(PROTEIN == "Q9QUI0")

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

#summming intensities 
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

psm.count.table <-psm.count.table%>%
  subset(Var1 %in% rownames(data_raw))
data_raw<-as.data.frame(data_raw)
#SL normalization

exp1_raw <- data_raw[c(1:4)]
exp2_raw <- data_raw[c(5:8)]
target <- mean(c(colSums(exp1_raw), colSums(exp2_raw)))
norm_facs <- target / colSums(exp1_raw)
exp1_sl <- sweep(exp1_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp2_raw)
exp2_sl <- sweep(exp2_raw, 2, norm_facs, FUN = "*")

data_sl <- cbind(exp1_sl, exp2_sl)

# Combat batch correction

design <- meta_group[-c(9:10),]
design$SampleGroup <- factor(design$SampleGroup, levels = c("Cul3_HET", "WT"))
mod <- model.matrix(~ SampleGroup + Gender, data = design)
batch <- design$Batch
data_sl=as.matrix(data_sl)
data_combat <- ComBat(dat = data_sl, batch = batch, mod = mod, par.prior = TRUE)
data_combat <-  data_combat[apply(data_combat, 1, function(x) all(x > 0)), ]
data_combat_A<-as.data.frame(data_combat)

#t-test

pval.372 = apply(log2(data_combat_A), 1, function(x) 
  t.test(as.numeric(x[c(1,2,3,4)]), as.numeric(x[c(5,6,7,8)]))$p.value)


logFC.372 = rowMeans(log2(data_combat_A[,c(1,2,3,4)]))-rowMeans(log2(data_combat_A[,c(5,6,7,8)]))
ttest.results = data.frame(gene=rownames(data_combat_A),
                           logFC=logFC.372,P.Value = pval.372, 
                           adj.pval = p.adjust(pval.372,method = "BH"))
write.csv(ttest.results, "ttest.results.csv",)

psm.count.table_ttest<-psm.count.table%>%
  filter(Var1 %in% ttest.results$gene)

ttest.results$PSMcount = psm.count.table_ttest[ttest.results$gene,]$Freq
ttest.results = ttest.results[with(ttest.results, order(adj.pval)), ]
head(ttest.results)
ttest.results %>% 
  filter(adj.pval <= 0.1 & abs(logFC) >= log2(1.5)) %>%
  count(adj.pval)

#ANOVA

mod_mat <- mod
colnames(mod_mat) <- make.names(colnames(mod_mat))

lm_fit <- lmFit(log2(data_combat), design = mod_mat)
contr_mat <- makeContrasts(
  HET = -SampleGroupWT,
  levels = colnames(mod_mat)
)

contr_fit <- contrasts.fit(lm_fit, contr_mat)

ord.t = contr_fit$coefficients[, 1]/contr_fit$sigma/contr_fit$stdev.unscaled[, 1]
ord.p = 2*pt(-abs(ord.t), contr_fit$df.residual)
ord.q = p.adjust(ord.p,method = "BH")
anova.results = data.frame(gene=names(contr_fit$sigma),
                           logFC=contr_fit$coefficients[,1],
                           t=ord.t, 
                           P.Value=ord.p, 
                           adj.P.Val = ord.q)

psm.count.table_anova<-psm.count.table%>%
  filter(Var1 %in% anova.results$gene)

anova.results$PSMcount = psm.count.table_anova[anova.results$gene,]$Freq
anova.results = anova.results[with(anova.results,order(P.Value)),]

head(anova.results)

anova.results %>% 
  filter(adj.P.Val <= 0.1 & abs(logFC) >= log2(1.5)) %>%
  count(gene)

#Limma

ebayes <- eBayes(contr_fit)
limma.results = topTable(ebayes,coef = "HET",n= Inf)
limma.results$gene = rownames(limma.results)
dpe <- limma.results %>%
  subset(adj.P.Val <= 0.1)

psm.count.table_limma<-psm.count.table%>%
  filter(Var1 %in% limma.results$gene)

limma.results$PSMcount = psm.count.table_limma[as.factor(limma.results$gene),]$Freq

head(limma.results)


dat.temp = data.frame(var = ebayes$sigma^2,
                      PSMcount = psm.count.table_limma[as.factor(names(ebayes$sigma)),]$Freq)

dat.temp.filter = dat.temp[dat.temp$PSMcount<40,]
op <- par(mfrow=c(1,2), mar=c(4,4,4,1), oma=c(0.5,0.5,0.5,0))
plot(log2(psm.count.table[rownames(data_combat_A),2]),
     log2(data_combat_A[,1])-log2(data_combat_A[,2]),pch=20,cex=0.5,
     xlab="log2(PSM count)",ylab="log2 ratio between two HET replicates")    

boxplot(log2(var)~PSMcount,dat.temp.filter,xlab="PSMcount",
        ylab="log2(Variance)")

ebayes$count = psm.count.table_limma[as.factor(rownames(ebayes$coefficients)),]$Freq


fit3 = spectraCounteBayes(ebayes)
DEqMS.results = outputResult(fit3,coef_col = 1)
head(DEqMS.results)
dpe_sp_ct <- DEqMS.results %>%
  subset(sca.adj.pval <= 0.1)

head(mouse_anno)
mouse_anno_sort_DEqMS=mouse_anno[match(rownames(DEqMS.results), mouse_anno$Accessions),]
DEqMS.results_annot<-cbind(DEqMS.results, mouse_anno_sort_DEqMS)
head(DEqMS.results_annot)
write.csv(DEqMS.results_annot, "DEqMS.results_annot_E17.5_CX.csv",)

## Comparison prior variance
plotFitCurve(fit3,type = "boxplot",xlab="PSM count")  
x = fit3$count
y = fit3$sigma^2

limma.prior = fit3$s2.prior
DEqMS.prior = fit3$sca.priorvar

plot(log2(x),log(y),pch=20, cex=0.5,xlab = "log2(PSMcount)",
     ylab="log(Variance)")
abline(h = log(limma.prior),col="green",lwd=3 )
k=order(x)
lines(log2(x[k]),log(DEqMS.prior[k]),col="red",lwd=3 )
legend("topright",legend=c("Limma prior variance","DEqMS prior variance"),
       col=c("green","red"),lwd=3)   

## Comparison posterior variance

x = fit3$count
y = fit3$s2.post

op <- par(mfrow=c(1,2), mar=c(4,4,4,1), oma=c(0.5,0.5,0.5,0))
plot(log2(x),log(y),pch=20, cex=0.5, xlab = "log2(PSMcount)", 
     ylab="log(Variance)",
     main="Posterior Variance in Limma")

y = fit3$sca.postvar
plot(log2(x),log(y),pch=20,cex=0.5, xlab = "log2(PSMcount)",
     ylab="log(Variance)", 
     main="Posterior Variance in DEqMS")

# P value distributions top 500 and 1000 genes

plot(sort(-log10(limma.results$P.Value),decreasing = TRUE)[1:500], 
     type="l",lty=2,lwd=2, ylab="-log10(p-value)",
     xlab="Proteins ranked by p-values",
     col="purple")
lines(sort(-log10(DEqMS.results$sca.P.Value),decreasing = TRUE)[1:500], 
      lty=1,lwd=2,col="red")
lines(sort(-log10(anova.results$P.Value),decreasing = TRUE)[1:500], 
      lty=2,lwd=2,col="blue")
lines(sort(-log10(ttest.results$P.Value),decreasing = TRUE)[1:500], 
      lty=2,lwd=2,col="orange")
legend("topright",legend = c("Limma","DEqMS","Anova","t.test"),
       col = c("purple","red","blue","orange"),lty=c(2,1,2,2),lwd=2)

#

library(ggplot2)
dat_filt_log<-dat_filt
dat_filt_log[, 4:13]<-log2(dat_filt[, 4:13])
dat.psm.log.ordered = dat_filt_log[,c("SEQUENCE","PROTEIN",
                                     sort(colnames(dat_filt_log)[4:11]))]
peptideProfilePlot(dat=dat.psm.log.ordered,col=2,gene="Q9JLV5")
Q9QUI0
peptideProfilePlot(dat=dat.psm.log.ordered,col=2,gene="Q9QUI0")
#correlation heatmaps

library( pheatmap )
cm <- cor( (data_combat_A) )
# rearrange columns so that same sample types are together
cm.order = cm[order(colnames(cm)),order(colnames(cm))]

pheatmap(cm.order,
         cluster_rows=T, cluster_cols=T,
         color = colorRampPalette(c("blue", "white", "red"))(100)
         )
dm <- as.matrix( dist( t( data_combat_A ) ))
dm.order = dm[order(colnames(cm)),order(colnames(cm))]
pheatmap(dm.order,
         cluster_rows=T, cluster_cols=T,clustering_method = "ward.D2",
         color = colorRampPalette(c("red", "white"))(100))  
