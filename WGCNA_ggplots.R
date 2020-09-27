setwd("/Users/pramod/Desktop/Clu3/WGCNA_Cul3/E17.5_WGCNA_CX/Prot_CX_E17.5_bicor_signed_sft18/")
library(tidyverse)
library(limma)
library(edgeR)
library(WGCNA)
library(ggplot2)
library(WGCNA)
library(ggplot2)
library(igraph)
# BiocManager::install("pSI")
library(pSI)
library(nlme)
library(gridExtra)
enableWGCNAThreads()
allowWGCNAThreads()
load("CX_norm_counts_bicor100_DeepSplit_0sft_power18.RData")
geneTree = networks$tree
datExpr=networks$datExpr
merged = networks$merged
modules = merged$colors
genes=rownames(datExpr)
MEs = networks$MEs
kMEtable = networks$kMEtable
CX_meta=readxl::read_xlsx("meta_group_wgcna.xlsx", sheet = 1)
datMeta=CX_meta[-c(9:10),]
datMeta$SampleGroup=factor(datMeta$SampleGroup,levels = c("WT","Cul3_HET"))
#datMeta=datMeta[match(colnames(datExpr),rownames(datMeta)),]
pSI.zhang = read.csv(file="Zhang_Human_Neuro2015_pSI_GSE21653.csv",row.names=1,stringsAsFactors = F)

pSI.zeisel = read.csv(file="Zeisel_level1_Mouse_Science2015_pSI.ensg.csv",row.names=1,stringsAsFactors = F)

pSI.goldmann = read.csv(file="Goldman_levelHybrid_Mouse_NatImmunol2015_pSI.ensg.csv",row.names=1,stringsAsFactors = F)

markers.Habib = read.csv("Habib_NatMeth2017_Human_DroncSeq_clusters.csv",row.names=1)
modTrait=read.delim("modTrait_0_100.txt")
library(plyr)
sft=18

for(i in 2:ncol(MEs$eigengenes)) {

  me = MEs$eigengenes[,i]
  moduleNumber = as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
  moduleColor = labels2colors(moduleNumber)
  moduleGenes = genes[modules==moduleNumber]
  
  pdf(file=paste0( moduleNumber, moduleColor, ".pdf"), width=20, height=11, useDingbats = FALSE)
  s = summary(lm(me ~ SampleGroup,data=datMeta))$coefficients
  dat2=data.frame()
  for(grp in c("Cul3_HET")) {
    rowID = paste0("SampleGroup", grp)
    dat2 = rbind.fill(dat2,
                 data.frame(Module=moduleColor, moduleNumber= moduleNumber, Group=grp,
                            beta = s[rowID, "Estimate"], SE = s[rowID, "Std. Error"], t=s[rowID, "t value"], p=s[rowID, "Pr(>|t|)"]))
  }
  dat2$fdr = modTrait$fdr[modTrait$moduleNumber == moduleNumber]
  dat2$p.symbol = ""; dat2$p.symbol[dat2$fdr<.05] ="*"
  g1=ggplot(dat2, aes(x=Group, y=beta, label=p.symbol)) + geom_bar(stat="identity",fill=moduleColor) +  
    geom_errorbar(aes(ymin=(beta - SE), ymax=(beta + SE)), width=0.25,size=0.25) +
    geom_text(color="red",size=8,aes(y=beta+ 1.3*sign(beta)*SE ))+
    ggtitle("Module-Trait Association")+
    xlab("")+
    ylab("Linear Regression Beta")
  
  dat=data.frame(Eigengene=me, Gender=datMeta$Gender,Group=datMeta$SampleGroup)
  g2=ggplot(dat,aes(x=Group,y=Eigengene,color=Group)) + geom_point(size=5)+    
    scale_color_manual(values = c("black","red","blue"))
  
  go=readxl::read_xlsx("GOEnrich_KME_sft18_signed_norm_counts_bicor_sig_moderate_0_MM_100.xlsx", sheet = paste0(i-1," ", moduleColor))
  go=go%>%
    top_n(15)
  g3=ggplot(go, aes(x=reorder(term.name, -log10(p.value)), y=-log10(p.value))) + geom_bar(stat="identity", fill="royalblue") + 
    coord_flip() + geom_hline(yintercept=-log10(0.05), lty=2, color="red")
  
  hubGenes = moduleGenes[order(kMEtable[moduleGenes,i], decreasing = T)[1:min(20,length(moduleGenes))]]
  annotation<-readxl::read_xlsx("E17.5_CX_combat_norm_data_log.xlsx", sheet = 1)
  annotation<-as.data.frame(annotation)
  hubGene.symbols=annotation[match(hubGenes, annotation$Protein_ID),]$Gene
  adjMat = adjacency(t(datExpr[hubGenes,]),type = "signed",corFnc = "cor", power=sft)
  adjMat[adjMat < quantile(adjMat,0.1)]=0
  graph <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=T,diag=F)
  plotcord= data.frame(layout_with_fr(graph))
  colnames(plotcord) = c("X1","X2")
  edgelist <- get.edgelist(graph,names = F)
  edges <- data.frame(plotcord[edgelist[,1],], plotcord[edgelist[,2],])
  colnames(edges) <- c("X1","Y1","X2","Y2")
  plotcord = cbind(plotcord, data.frame(gene=hubGene.symbols))
  
  g4=ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = 0.5, colour="grey") + 
    geom_point(aes(X1, X2), data=plotcord,color=moduleColor,size=4) + geom_text(aes(x=X1, y=X2+.2, label=gene),fontface="bold",size=6,data=plotcord) +
    theme_classic()
  
  
  grid.arrange(grobs=list(g1,g2,g3,g4), layout_matrix=rbind(c(1,1,2,2,2,5,5),c(3,3,3,4,4,4,4)),
               top=paste0("Module ", moduleNumber, " (", moduleColor, ")"),padding=unit(2, "line"))
               dev.off()
}
  


