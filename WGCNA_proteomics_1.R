setwd("/Users/pramod/Desktop/Clu3/WGCNA_Cul3/E17.5_WGCNA_CX/Prot_CX_E17.5_bicor_signed_sft18/")
library(tidyverse)
library(limma)
library(edgeR)
library(WGCNA)
enableWGCNAThreads()
allowWGCNAThreads()
TPM_CX<-readxl::read_xlsx("E17.5_CX_combat_norm_data_log.xlsx", sheet = 1)
CX_meta=readxl::read_xlsx("meta_group_wgcna.xlsx", sheet = 1)
TPM_CX<-as.data.frame(TPM_CX)
rownames(TPM_CX)<-TPM_CX$Protein_ID; TPM_CX$Gene=NULL; TPM_CX$Protein_ID=NULL
# rownames(TPM_CX)=TPM_CX$Gene
# TPM_CX=TPM_CX[,match(rownames(CX_meta),colnames(TPM_CX))]
# design=model.matrix(~Genotype + Sex + W_1 + W_2 + W_3 + W_4 + W_5+W_6+W_7+W_8,
#                     data = CX_meta)
# dge.voom = voom(calcNormFactors(DGEList(counts = TPM_CX),
#                                 method = 'TMM'), design, plot = T)
# Y = as.matrix(dge.voom$E)
# X = as.matrix(design)
# beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
# regressed_data = Y - t(X[,c(3:ncol(X))] %*% beta[c(3:nrow(beta)),])

regressed_data<-TPM_CX

pdf("PCA_CX_regressed_counts.pdf", height = 8, width =8 )
pr<-prcomp(t(regressed_data[, -9]), center = TRUE, scale = TRUE)
pr$x %>%
  as.data.frame()%>%
  mutate(SampleName = rownames(.))%>%
  inner_join(CX_meta, by = "SampleName")%>%
  ggplot(
    data=.,
    aes(PC1, PC2, label = Sample, color = SampleGroup
    )
  )+
  geom_point(size = 4)+
  geom_text(vjust = 1)
dev.off()


powers = c(seq(1,9,by=1),seq(10,30,by=2)) 
sft_cor = pickSoftThreshold(data= t(regressed_data), 
                              networkType = "signed", corFnc="bicor",verbose=2,
                              powerVector=powers,blockSize = 20000)


pdf("sft_cx_bicor_regressed_counts_signed.pdf", height = 8, width = 8)
par(mfrow=c(1,2))
plot(sft_cor$fitIndices[,1], 
     -sign(sft_cor$fitIndices[,3])*sft_cor$fitIndices[,2],
     xlab="Soft Thresh Power", ylab="Scale free R^2",type="n")
text(sft_cor$fitIndices[,1], 
     -sign(sft_cor$fitIndices[,3])*sft_cor$fitIndices[,2], 
     labels = powers, cex = 0.7, col="red",  xlab="Soft Thresh Power", 
     ylab="Scale free R^2")
abline(h=0.8, col="black")
plot(sft_cor$fitIndices[,1], sft_cor$fitIndices[,5], 
     xlab = "Soft threshold power", ylab = "Mean connectivity", type = "n")
text(sft_cor$fitIndices[,1], sft_cor$fitIndices[,5], 
     labels = powers, cex = 0.7, col="black")
dev.off()

pdf("sft_bicor_18_signed.pdf", height = 8, width = 8)
net = blockwiseModules(datExpr=t(regressed_data[, -9]), maxBlockSize=20000,networkType="signed",corType="bicor",  
                       power = 18, mergeCutHeight= 0.1, minModuleSize= 40, pamStage=TRUE, reassignThreshold=1e-6, 
                       saveTOMFileBase="CX_WGCNA_signed_bicor_norm_counts_18", saveTOMs=TRUE, 
                       verbose = Inf, deepSplit=2)
dev.off()
save(file = "cx_e175_w8_net_bicor_norm_counts_ds_2_sft_18.RData",net)
load("./CX_WGCNA_signed_bicor_norm_counts_18-block.1.RData")
sft=18
for (mm in c(20,30,40,50,70,100, 150, 200, 250, 300)){
  print(mm)
  for (ds in seq(0,4)){
    print(ds)
    ds = ds; minModSize = mm; dthresh = 0.1; pam = FALSE
    networks=list()
    networks$datExpr=regressed_data
    networks$tree = hclust(1-TOM, method="average")
    networks$cut = cutreeHybrid(dendro = networks$tree, pamStage=pam, minClusterSize= minModSize, 
                                cutHeight = 0.99999, deepSplit=ds, distM=as.matrix(1-TOM))
    networks$merged= mergeCloseModules(exprData= t(networks$datExpr), colors = networks$cut$labels, cutHeight=dthresh)
    networks$MEs = moduleEigengenes(t(networks$datExpr), colors=networks$merged$colors, softPower=sft)
    networks$kMEtable = signedKME(t(networks$datExpr), datME = networks$MEs$eigengenes,corFnc = "bicor")
    tryCatch(
      {
        pdf(paste("./CX_norm_counts_bicor_",mm, "_DeepSplit_",ds,"sft_power",sft,".pdf",sep = ""), height = 8, width = 8)
        plotDendroAndColors(networks$tree, colors=labels2colors(networks$merged$colors), dendroLabels = F)
        dev.off()
        
      },
      error = function(e) {
        message("Too many modules to plot, skipping...")
      }
    )
    
    save(file=paste("./CX_norm_counts_bicor",mm, "_DeepSplit_",ds,"sft_power",sft,".RData",sep = ""), networks)
  }
}

CX_meta<-CX_meta[-c(9:10),]
library(plyr)
# module-trait plots
pdf("CX_plots_module_trait_counts_norm_sft_18_bicor_sig_A.pdf")
for (mm in c(20,30,40,50,70,100,150,200, 250, 300)){
  print(mm)
  for (ds in seq(0,4)){ 
    print(ds)
    load(paste("./CX_norm_counts_bicor",mm, "_DeepSplit_",ds,"sft_power",sft,".RData",sep = ""))
    datExpr=networks$datExpr
    geneTree=networks$tree
    merged=networks$merged
    modules=merged$colors
    print(length(unique(modules))-1)
    MEs=networks$MEs
    kMEtable=networks$kMEtable
    table(rownames(CX_meta) == colnames(regressed_data))
    CX_meta$SampleGroup=factor(CX_meta$SampleGroup,levels = c("WT","Cul3_HET"))
    
    modTrait=data.frame()
    for(i in 2:length(unique(modules))) {
      me = MEs$eigengenes[,i]
      moduleNumber = as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
      moduleColor = labels2colors(moduleNumber)
      s = summary(lm(me ~ SampleGroup + Gender,data=CX_meta))$coefficients
      for(grp in c("Cul3_HET")){
        
        rowID = paste0("SampleGroup", grp)
        
        modTrait = rbind.fill(modTrait,
                              data.frame(Module=moduleColor, moduleNumber= moduleNumber, Group=grp,
                                         beta = s[rowID, "Estimate"], SE = s[rowID, "Std. Error"], t=s[rowID, "t value"], p=s[rowID, "Pr(>|t|)"]))}
      
      modTrait$fdr=p.adjust(modTrait$p,method = "fdr")
      modTrait$signedLog10fdr = -log10(modTrait$fdr) * sign(modTrait$beta)
      modTrait$signedLog10fdr[modTrait$fdr > .05] = 0
      modTrait$text = signif(modTrait$beta, 1)
      print(modTrait)
      write.table(modTrait, file="modTrait_0_100.txt", sep = "\t",)
      modTrait_sig<-modTrait %>%
        filter(!signedLog10fdr == 0)
      print(modTrait_sig)
      # modTrait$text[modTrait$fdr > 0.05] = ""
      # print(ggplot(modTrait, aes(x=Module,y=Group, label=text)) +
      #         geom_tile(aes(fill=signedLog10fdr),color="grey60") +
      #         scale_fill_gradient2(low = "blue", high = "red","[beta]\nsigned\n-log10FDR\n") +
      #         geom_text(size=3, color="black")  +
      #         ggtitle(paste("soft_power_", sft,"_module_trait_module_size_", mm, "_DeepSplit_",ds, sep="")))
      
      
      # print(ggplot(modTrait,aes(x=Module,y=beta, fill=Group)) + geom_bar(stat="identity", position=position_dodge(width=.75),width=.75) +
      #   geom_errorbar(aes(x=Module,ymin=beta-SE, ymax=beta+SE),width=.25,position=position_dodge(width = .75)) + 
      #   theme(axis.text.x = element_text(angle=-60, hjust=0)) + labs(y="log2FC", x=""))
    }
    print(ggplot(modTrait_sig, aes(x=Module,y=Group, label=text)) +
            geom_tile(aes(fill=signedLog10fdr),color="grey60") +
            scale_fill_gradient2(low = "blue", high = "red","[beta]\nsigned\n-log10FDR\n") +
            geom_text(size=3, color="black")  +
            ggtitle(paste("soft_power_", sft,"_module_trait_module_size_", mm, "_DeepSplit_",ds, sep="")))
    
  }}
dev.off()

library(gProfileR)
for (mm in c(20,30, 40, 50, 70, 100, 150, 200, 250, 300)){
  print(mm)
  for (ds in seq(0,4)){
    print(ds)
    load(paste("./CX_norm_counts_bicor",mm, "_DeepSplit_",ds,"sft_power",sft,".RData",sep = ""))
    datExpr=networks$datExpr
    genes=rownames(datExpr)
    merged=networks$merged
    modules=merged$colors
    MEs = networks$MEs
    kMEtable=networks$kMEtable
    
    unique_module_labels <- sort(unique(modules))
    unique_module_labels <- unique_module_labels[unique_module_labels != 0]
    unique_module_colours <- WGCNA::labels2colors(unique_module_labels)
    
    GO_enrich<-list()
    for (i in seq_along(unique_module_labels)) {
      for (filter in c("strong", "none", "moderate")){
      moduleNumber = unique_module_labels[[i]]
      moduleColor = unique_module_colours[[i]]
      moduleGenes = genes[modules==moduleNumber]
      #moduleGenes=moduleGenes[order(kMEtable[moduleGenes,i], decreasing = T)],
      annotation<-readxl::read_xlsx("E17.5_CX_combat_norm_data_log.xlsx", sheet = 1)
      annotation<-as.data.frame(annotation)
      moduleGenes=annotation[match(moduleGenes, annotation$Protein_ID),]$Gene
      go = gprofiler(query=moduleGenes, 
                     organism = "mmusculus",
                     correction_method = "fdr",
                     exclude_iea = T,
                     hier_filtering = filter,
                     custom_bg = genes, 
                     src_filter = c("GO:BP","GO:MF"),
                     ordered_query = F)%>%
      arrange(p.value)%>%
      filter(overlap.size > 1)
      module_name <- paste(moduleNumber, moduleColor)
      GO_enrich[[module_name]] <- go
    }
    writexl::write_xlsx(
      GO_enrich, 
      path = paste("GOEnrich_sft18_signed_norm_counts_bicor_",filter,"_",ds,"_MM_",mm,".xlsx", sep = "")
      
    )
    
    }
  }
}
