
  # 0. loading library
  
  ```{r}
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(clustree)
library(infercnv)
library(phylogram)
library(ape)
library(org.Hs.eg.db)
library(RColorBrewer)
library(SpaCET)
```

# quality control, filter, merge and cluster -------------------------------------------------------

new.C15 <- PercentageFeatureSet(C15, pattern = "^MT-",col.name = "percent.mt")
new.C16 <- PercentageFeatureSet(C16, pattern = "^MT-",col.name = "percent.mt")
new.C17 <- PercentageFeatureSet(C17, pattern = "^MT-",col.name = "percent.mt")


# Exclude mitochondrial genes from the counts matrix

filtered_C15 <- new.C15[!grepl("^MT-", rownames(new.C15)), ]
filtered_C15 <- filtered_C15[!grepl("^RP[SL]", rownames(filtered_C15)), ]
filtered_C15 <- subset(filtered_C15,subset = nCount_Spatial >300)
filtered_C15 <- subset(filtered_C15,subset = nFeature_Spatial >1000)
filtered_C16 <- new.C16[!grepl("^MT-", rownames(new.C16)), ]
filtered_C16 <- filtered_C16[!grepl("^RP[SL]", rownames(filtered_C16)), ]
filtered_C16 <- subset(filtered_C16,subset = nCount_Spatial >300)
filtered_C16 <- subset(filtered_C16,subset = nFeature_Spatial >1000)
filtered_C17 <- new.C17[!grepl("^MT-", rownames(new.C17)), ]
filtered_C17 <- filtered_C17[!grepl("^RP[SL]", rownames(filtered_C17)), ]
filtered_C17 <- subset(filtered_C17,subset = nCount_Spatial >300)
filtered_C17 <- subset(filtered_C17,subset = nFeature_Spatial >1000)

table(filtered_C15$label2)
table(filtered_C16$label2)
table(filtered_C17$label2)


sct15 <- SCTransform(filtered_C15, assay = "Spatial",return.only.var.genes = T,do.scale = T,vars.to.regress = "percent.mt" ,verbose = FALSE) 
sct16 <- SCTransform(filtered_C16, assay = "Spatial",return.only.var.genes = T,do.scale = T,vars.to.regress = "percent.mt" ,verbose = FALSE) 
sct17 <- SCTransform(filtered_C17, assay = "Spatial",return.only.var.genes = T,do.scale = T,vars.to.regress = "percent.mt" ,verbose = FALSE) 

```

# 3. merge and find cluster


### no harmony no batch corrected

sct.3 <- merge(sct15,merge(sct16,sct17))
VariableFeatures(sct.3) <- c(VariableFeatures(sct15),VariableFeatures(sct16),VariableFeatures(sct17)) 
table(duplicated(VariableFeatures(sct.3)))  
sct.3 <- RunPCA(sct.3,assay = "SCT")
DimPlot(sct.3,reduction = "pca") # pca clustered 

rm(sct15,sct16,sct17,filtered_C15,filtered_C16,filtered_C17)

sct.3  <- FindNeighbors(sct.3,assay = "SCT",dims = 1:30)
sct.3 <-  FindClusters(sct.3,resolution= 1.2)
sct.3 <- RunUMAP(sct.3,dims = 1:30,assay = "SCT")

DimPlot(sct.3,reduction = "umap",label = T,label.box = T,pt.size = 1)
DimPlot(sct.3,reduction = "umap",group.by = "label2",label = T,label.box = T,pt.size = 1)
DimPlot(sct.3,reduction = "umap",group.by = "orig.ident",label = T,label.box = T,pt.size = 1)
SpatialDimPlot(sct.3,label = T)
SpatialDimPlot(sct.3,group.by = "label2",label = T)


# # find all markers  -------------------------------------------------------
# 
DefaultAssay(sct.3) <- "Spatial"
sct.3 <- SetIdent(sct.3,value = "seurat_clusters")

sct.3.allmarkers <- FindAllMarkers(sct.3,only.pos = FALSE, min.pct = 0.5, logfc.threshold = 0.5)
sct.3.allmarkers$p_val <- round(sct.3.allmarkers$p_val,5)
sct.3.allmarkers$p_val_adj <- round(sct.3.allmarkers$p_val_adj,10)
sct.3.sigmarkers <- sct.3.allmarkers %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_log2FC)

sct.3.allmarkers$diff <- sct.3.allmarkers$pct.1-sct.3.allmarkers$pct.2
sct.3.diffmarkers <- sct.3.allmarkers %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = diff)

sct.3.mergemarker

SpatialDimPlot(sct.3)

sct.3 <- NormalizeData(sct.3,assay = "Spatial",normalization.method = "LogNormalize")
sct.3 <- ScaleData(sct.3,features = rownames(sct.3))

DoHeatmap(sct.3, features = sct.3.diffmarkers$gene) + NoLegend() ## landscape 15 * 20
DoHeatmap(sct.3, features = sct.3.sigmarkers$gene) + NoLegend()  ## landscape 15 * 20

# 

```

# 4.define for the cluster


new.cluster.ids <- c("PDAC", 
                     "Acinar", "Stroma", "Acinar", "Stroma", "Stroma",
                     "Stroma","Immune","IPMN","Stroma","Ductal",
                     "Stroma","Ductal","Immune","Islet","Immune")
names(new.cluster.ids) <- levels(sct.3)
new_sct.3 <- RenameIdents(sct.3, new.cluster.ids)
new_sct.3$new.cluster.ids <- new_sct.3@active.ident
Idents(new_sct.3) <- "new.cluster.ids"



# Generate a discrete color palette with 8 colors
new.cluster.ids_color <- brewer.pal(7, "Set2")
names(new.cluster.ids_color) <- levels(new_sct.3)
# Plot the colors
par(mar = c(3, 3, 2, 2))
plot(1:7, rep(1, 7), pch = 19, cex = 4, col = new.cluster.ids_color, xlab = "", ylab = "", xaxt = 'n', yaxt = 'n', bty = 'n')+
  text(1:8, rep(1, 8) - 0.05, labels = new.cluster.ids_color, pos = 1, cex = 0.8, col = "black")+
  text(1:8, rep(1, 8) - 0.15, labels = names(new.cluster.ids_color), pos = 1, cex = 0.8, col = "black")


DimPlot(new_sct.3, reduction = "umap", label = TRUE, pt.size = 0.5,cols = new.cluster.ids_color)
DimPlot(new_sct.3,group.by = "label2",label = T,pt.size = 0.5)
DimPlot(new_sct.3,group.by = "orig.ident",label = T,pt.size = 0.5)

DimPlot(new_sct.3,split.by = "orig.ident",cols = new.cluster.ids_color,label = T,pt.size = 0.5)

SpatialDimPlot(new_sct.3,images = "inmage",crop = F, pt.size.factor = 7000,cols = new.cluster.ids_color)
SpatialDimPlot(new_sct.3,images = "inmage.1",crop = F, pt.size.factor = 7000,cols = new.cluster.ids_color)
SpatialDimPlot(new_sct.3,images = "inmage.1.1",crop = F, pt.size.factor = 7000,cols = new.cluster.ids_color)

SpatialDimPlot(new_sct.3,crop = F, pt.size.factor = 7000,group.by = "label2")


SpatialFeaturePlot(new_sct.3,images = "inmage.1.1",features = "TSPAN8",crop = F, pt.size.factor = 7000)

SpatialDimPlot(new_sct.3,group.by = "label2",label = T,crop = F, pt.size.factor = 7000)
SpatialDimPlot(new_sct.3,group.by = "orig.ident",label = T,crop = F, pt.size.factor = 7000)


```

## 1. heatmap or vlnplot
```{r}
new_sct.3.markers <- FindAllMarkers(new_sct.3, only.pos = TRUE)

new_sct.3.markers$diff <- new_sct.3.markers$pct.1 - new_sct.3.markers$pct.2
new_sct.3.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1 & pct.1 >0.9)%>%
  slice_head(n = 20) %>%
  ungroup() -> top20


new_sct.3.markers %>%
  group_by(cluster) %>%
  dplyr::filter(diff > 0.5 & pct.1 >0.9)%>%
  slice_head(n = 20) %>%
  ungroup() -> top20_diff


merge_markers <- merge(top20[,c(2,6,7)],top20_diff[,c(6,7,8)],by = "gene")
DotPlot(new_sct.3,features = unique(top20$gene))+coord_flip()+ NoLegend()
DotPlot(new_sct.3,features = unique(top20_diff$gene))+coord_flip()+ NoLegend()
DotPlot(new_sct.3,features = c("MUC5AC","AGR3",
                               "AMY2A","CELA3A",
                               "VIM","COL1A1",
                               "IGKC","IGHG4",
                               "TFF2","MUC1",
                               "CFTR","PLA2G1B",
                               "ABCC8","IAPP"))+coord_flip()+ NoLegend()



for (i in c("MUC5AC","AGR3","MUC6","MUC2",
            "AMY2A","CELA3A",
            "VIM","COL1A1",
            "IGKC","IGHG4",
            "TFF2","MUC1",
            "CFTR","KRT19",
            "INS","IAPP")) {
  p1 = SpatialFeaturePlot(new_sct.3,features = i,crop = F, pt.size.factor = 7000)
  print(p1)
  
}

VlnPlot(new_sct.3,features = c("MUC5AC","AGR3",
                               "AMY2A","CELA3A",
                               "VIM","COL1A1",
                               "IGKC","IGHG4",
                               "TFF2","MUC1",
                               "CFTR","KRT19",
                               "INS","IAPP"),
        cols = new.cluster.ids_color,stack = T,flip = T,assay = "Spatial",slot = "data")


## 2. spatial genes plot
for (i in c("MUC5AC","AGR3","MUC6","MUC2",
            "AMY2A","CELA3A",
            "VIM","COL1A1",
            "IGKC","IGHG4",
            "TFF2","MUC1",
            "CFTR","KRT19",
            "INS","IAPP")) {
  p1 = SpatialFeaturePlot(new_sct.3,features = i,crop = F, pt.size.factor = 7000)
  print(p1)
  
}

## 3. vlnplot

a = subset(new_sct.3,subset = new.cluster.ids %in% c("Ductal","IPMN","PDAC"))
a$new.cluster.ids <- as.factor(as.character(a$new.cluster.ids))
a <- SetIdent(a, value = a$new.cluster.ids)

# 重新设置ident顺序
a@active.ident <- factor(a@active.ident, levels = c("Ductal", "IPMN", "PDAC"))

# 检查新的ident顺序
levels(a@active.ident)

for (i in c("TFF2","MUC6","PGC","AQP5","ONECUT3","SOX9","GATA6",
            "KRT17","AREG","AGR3","GATA3","ELF3","KLF4","HMGA1")) {
  p1 = VlnPlot(a,features = i,
               cols = c("PDAC"= "#00BFC4","IPMN"= "#F8766D","Ductal"= "#00BA38"))
  print(p1)
  
}



## 4.overlap of annotate vs cluster

a = as.data.frame(table(new_sct.3$new.cluster.ids,new_sct.3$label2))

proportion.a <- a %>% 
  group_by(Var2) %>%
  mutate(percentage = Freq / sum(Freq) * 100)

# Manually set the order of levels for Var2
proportion.a$Var2 <- factor(proportion.a$Var2, levels = c("Invasive_PDAC","Immune_rich_region","Stroma","HG_IPMN","Immune_infiltration",
                                                          "Pancreatic_duct","Normal_acinar","Pancreatic_islet"))
ggplot(proportion.a, aes(x = Var2, y = Var1)) +
  geom_point(aes(size = percentage, color = Var1)) +
  scale_size(range = c(0, 10)) +
  scale_color_manual(values = new.cluster.ids_color) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Annotation")+
  ylab("Cluster")


## 5.proportion of annotate/cluster sample

a = as.data.frame(table(new_sct.3$orig.ident,new_sct.3$label2))


a = as.data.frame(table(new_sct.3$orig.ident,new_sct.3$new.cluster.ids))

proportion.a <- a %>% 
  group_by(Var1) %>%
  mutate(percentage = Freq / sum(Freq) * 100)

print(proportion.a)
write.csv(proportion.a,file = "proportion_sample_cluster.csv")
pdf("rm_mt_rp_genes_results/BarPlot_sample_NewClusterID_overlap.pdf",height = 8,width = 5)
ggplot(proportion.a, aes(fill = Var2, y = percentage, x = Var1)) +
  geom_bar(stat = "identity") +  
  scale_fill_manual(values = new.cluster.ids_color) +
  ylab("Percentage of Cluster") +  
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1))+
  theme_classic()
dev.off()

rm(a,proportion.a)
```

# 5.diff analysis across Ductal \> IPMN \> PDAC

## 1.filter all overlap ductal/IPMN/PDAC

```{r}

ductal_tumor <- subset(new_sct.3,subset = new.cluster.ids %in% c("Ductal","IPMN","PDAC"))
# 1351 spots
table(ductal_tumor$label2,ductal_tumor$new.cluster.ids)
table(ductal_tumor$orig.ident,ductal_tumor$new.cluster.ids)

# second filter by pathologist label
ductal_tumor <- subset(ductal_tumor,subset = label2 %in% c("Pancreatic_duct","HG_IPMN","Invasive_PDAC"))
table(ductal_tumor$label2,ductal_tumor$new.cluster.ids)
table(ductal_tumor$orig.ident,ductal_tumor$new.cluster.ids)

```

Annotated IPMN has IPMN and Ductal Annotated PDAC has PDAC, IPMN, Ductal

In total: Ductal = 374 spots IPMN = 315 spots PDAC = 578 spots

## 2. find deg

```{r}
deg_Ductal <- FindMarkers(ductal_tumor,ident.1 = c("PDAC","IPMN"),ident.2 = "Ductal")
deg_Ductal$symbol <- rownames(deg_Ductal)
deg_Ductal <- deg_Ductal[order(deg_Ductal$avg_log2FC,decreasing = T),]
deg_PDAC_IPMN <- FindMarkers(ductal_tumor,ident.1 = "PDAC",ident.2 = "IPMN")
deg_PDAC_IPMN$symbol <- rownames(deg_PDAC_IPMN)
deg_PDAC_IPMN <- deg_PDAC_IPMN[order(deg_PDAC_IPMN$avg_log2FC,decreasing = T),]

merge_deg <- merge(deg_PDAC_IPMN,deg_Ductal,by = "symbol")
merge_deg <- filter(merge_deg,merge_deg$p_val_adj.x<0.05 & merge_deg$p_val_adj.y < 0.05 )

merge_deg <- merge_deg[,c(1,3,6,8,10)]




tf.list <- read.csv("./Genesets/Animal_v3_Homo_sapiens_TF.csv")

merge_tf <- merge_deg[merge_deg$symbol %in% tf.list$Symbol,]



write.csv(merge_deg,file = "merge_deg.csv")
write.csv(merge_tf,file = "merge_tf.csv")

```

### 2.1 plot



merge_deg$group <- ifelse(merge_deg$avg_log2FC.y>2,
                          ifelse(merge_deg$avg_log2FC.x>1, "PDAC",
                                 ifelse(merge_deg$avg_log2FC.x< -1,"IPMN",NA)),
                          ifelse(merge_deg$avg_log2FC.y < -2,"Ductal",NA))
merge_deg$label <- ifelse(merge_deg$group %in% c("PDAC","IPMN","Ductal"),merge_deg$symbol,NA)

merge_tf$group <- ifelse(merge_tf$avg_log2FC.y>0,
                         ifelse(merge_tf$avg_log2FC.x>0, "PDAC",
                                ifelse(merge_tf$avg_log2FC.x< 0,"IPMN",NA)),
                         ifelse(merge_tf$avg_log2FC.y < 0,"Ductal",NA))
merge_tf$label <- ifelse(merge_tf$group %in% c("PDAC","IPMN","Ductal"),merge_tf$symbol,NA)


pdf("rm_mt_rp_genes_results/PointPlot_deg_detf.pdf",width = 20,height = 15)
ggplot(merge_deg,aes(x = avg_log2FC.x,y = avg_log2FC.y,color = group,size = p_val_adj.x))+
  geom_point()+
  theme_classic()+
  geom_vline(xintercept = 0)+
  geom_vline(xintercept = 1,linetype = "dashed")+
  geom_vline(xintercept = -1,linetype = "dashed")+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = 2,linetype = "dashed")+
  geom_hline(yintercept = -2,linetype = "dashed")+
  annotate(geom = "text",x = 5,y = 7,label = "PDAC only")+
  annotate(geom = "text",x = -5,y = 7,label = "IPMN only")+
  annotate(geom = "text",x = -5,y = -7,label = "Ductal, IPMN > PDAC")+
  annotate(geom = "text",x = 5,y = -7,label = "Ductal, PDAC > IPMN")+
  xlab("avgLFC of PDAC vs IPMN")+
  ylab("avgLFC of Tumor vs Ductal")+
  xlim(-10,10)+
  geom_text_repel(data = subset(merge_deg, !is.na(label)), aes(label = label), size = 5,
                  point.padding = 0.5,
                  max.overlaps = 20)  # Exclude NA hits



ggplot(merge_tf,aes(x = avg_log2FC.x,y = avg_log2FC.y,color = group,size = p_val_adj.x))+
  geom_point()+
  theme_classic()+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  annotate(geom = "text",x = 4.5,y = 4.5,label = "PDAC only")+
  annotate(geom = "text",x = -4.5,y = 4.5,label = "IPMN only")+
  annotate(geom = "text",x = -4.5,y = -4.5,label = "Ductal, IPMN > PDAC")+
  annotate(geom = "text",x = 4.5,y = -4.5,label = "Ductal, PDAC > IPMN")+
  xlab("avgLFC of PDAC vs IPMN")+
  ylab("avgLFC of Tumor vs Ductal")+
  xlim(-5,5)+
  geom_text_repel(data = subset(merge_tf, !is.na(label)), aes(label = label), size = 5,
                  point.padding = 0.5,
                  max.overlaps = 20) 

dev.off()


pdf(paste0("rm_mt_rp_genes_results/sFeaturePlot_tf.pdf"),height = 5,width = 9)
for (i in merge_tf$label) {
  p1 = SpatialFeaturePlot(new_sct.3,features = i,crop = F, pt.size.factor = 7000)
  print(p1)
  
}
dev.off()



```

## 3.PDAC side enrich +/+
```{r}
library(clusterProfiler)
library(msigdbr)
library(openxlsx)
gmt_all_subtype <- read.xlsx("./Genesets/all_subtype.xlsx",sheet = "all_geneset")
gmt_go <- read.gmt("./Genesets/GSEA_GO_BP_CC_MF.gmt")
gmt_kegg <- read.gmt("./Genesets/GSEA_KEGG.gmt")
gmt_hm <- read.gmt("./Genesets/GSEA_Hallmark.gmt")


pdac_genes <- merge_deg[merge_deg$avg_log2FC.x>0 & merge_deg$avg_log2FC.y>0,]$symbol


## enrichR
library(enrichR)

dbs <- listEnrichrDbs()
dbs <- c("WikiPathway_2023_Human","MSigDB_Oncogenic_Signatures","Metabolomics_Workbench_Metabolites_2022","MSigDB_Hallmark_2020","KEGG_2021_Human","GO_Molecular_Function_2023","GO_Cellular_Component_2023","GO_Biological_Process_2023","Enrichr_Submissions_TF-Gene_Coocurrence")
enriched <- enrichr(pdac_genes, dbs)


pdf("rm_mt_rp_genes_results/Barplot_enrichr_PDACsite.pdf",width = 8,height = 8)
for (i in c(1:length(dbs))) {
  p1 = plotEnrich(enriched[[i]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",title = paste0(dbs[i]))
  print(p1)
}
dev.off()



```




## 4.IPMN side enrich +/-
```{r}
ipmn_genes <- merge_deg[merge_deg$avg_log2FC.x<0 & merge_deg$avg_log2FC.y>0,]$symbol

dbs <- c("WikiPathway_2023_Human","MSigDB_Oncogenic_Signatures","Metabolomics_Workbench_Metabolites_2022","MSigDB_Hallmark_2020","KEGG_2021_Human","GO_Molecular_Function_2023","GO_Cellular_Component_2023","GO_Biological_Process_2023","Enrichr_Submissions_TF-Gene_Coocurrence")
enriched <- enrichr(ipmn_genes, dbs)

pdf("rm_mt_rp_genes_results/Barplot_enrichr_IPMNsite.pdf",width = 8,height = 8)
for (i in c(1:length(dbs))) {
  p1 = plotEnrich(enriched[[i]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",title = paste0(dbs[i]))
  print(p1)
}
dev.off()


```


## 5. addmodule score

```{r}
library(clusterProfiler)
all_gmt <- openxlsx::read.xlsx("./Genesets/all_subtype.xlsx",sheet = "all_geneset")
hm.gmt <- read.gmt("./Genesets/GSEA_Hallmark.gmt")

new_sct.3 <- AddModuleScore(new_sct.3,
                            features = list(all_gmt[all_gmt$term == "Moffitt_Basal",2],
                                            all_gmt[all_gmt$term == "Moffitt_Classical",2]),
                            ctrl  = 100,
                            name = c("Moffitt_Basal_Score",
                                     "Moffitt_Classical_Score"))
for (i in unique(as.character(hm.gmt$term))) {
  new_sct.3 <- AddModuleScore(new_sct.3,
                              features = list(hm.gmt[hm.gmt$term == i,2]),
                              ctrl  = 100,
                              name = i)
}

pdf("./rm_mt_rp_genes_results/sFeaturePlot_SubtypeScore.pdf",width = 8,height = 8)

SpatialFeaturePlot(new_sct.3,features = colnames(new_sct.3@meta.data)[c(12:13)],keep.scale = "all",images = c("inmage","inmage.1"),crop = F,pt.size.factor = 7000)

dev.off()


pdf("./rm_mt_rp_genes_results/sFeaturePlot_Hallmark_scores.pdf",width = 9,height = 5)
for (i in colnames(new_sct.3@meta.data)[c(14:63)]) {
  p1 = SpatialFeaturePlot(new_sct.3,features = i,images = c("inmage","inmage.1"),keep.scale = "all",crop = F,pt.size.factor = 7000)
  print(p1)
}
dev.off()

pdf("./rm_mt_rp_genes_results/sFeaturePlot_Hallmark_rm_non_malignant_scores.pdf",width = 9,height = 5)
for (i in colnames(new_sct.3@meta.data)[c(12:61)]) {
  p1 = SpatialFeaturePlot(subset(new_sct.3,subset = new.cluster.ids %in% c("IPMN","PDAC")),features = i,images = c("inmage","inmage.1"),keep.scale = "all",crop = F,pt.size.factor = 7000)
  print(p1)
}
dev.off()

```
## 6. iCAF, mCAF score
```{r}
CAF.gmt <- read.csv("./Genesets/iCAF_myCAF_PMC6727976.csv")

new_sct.3 <- AddModuleScore(new_sct.3,
                            features = list(CAF.gmt[CAF.gmt$term == "iCAF",2],
                                            CAF.gmt[CAF.gmt$term == "myCAF",2]),
                            ctrl  = 100,
                            name = c("iCAF_Score",
                                     "myCAF_Score"))


pdf("./sFeaturePlot_CAF_Score.pdf",width = 8,height = 8)

SpatialFeaturePlot(new_sct.3,features = colnames(new_sct.3@meta.data)[c(12:13)],keep.scale = "all",images = c("inmage","inmage.1"),crop = F,pt.size.factor = 7000)

dev.off()

pdf("./VlnPlot_CAF_ImmuneMarkers_Score.pdf",width = 4,height = 3)

for (i in c(colnames(new_sct.3@meta.data)[c(12:13)],
            "CD40","CD4","CD8A","LYZ","CD14","CD19")) {
  p1 = VlnPlot(new_sct.3,features = i,idents = c("PDAC","IPMN"))
  print(p1)
}

dev.off()


pdf("./sFeaturePlot_ImmuneMarkers.pdf",width = 9,height = 5)
for (i in c("CD40","CD4","CD8A","LYZ","CD14","CD19")) {
  p1 = SpatialFeaturePlot(new_sct.3,features = i,images = c("inmage","inmage.1"),keep.scale = "all",crop = F,pt.size.factor = 7000)
  print(p1)
}
dev.off()

```

## 7. iCAF, mCAF score only stroma region
```{r}
CAF.gmt <- read.csv("./Genesets/iCAF_myCAF_PMC6727976.csv")

subset_new_sct.3 <- AddModuleScore(subset(new_sct.3,subset = new.cluster.ids %in% c("Stroma") & orig.ident %in% c("C15","C16")),
                                   features = list(CAF.gmt[CAF.gmt$term == "iCAF",2],
                                                   CAF.gmt[CAF.gmt$term == "myCAF",2]),
                                   ctrl  = 100,
                                   name = c("iCAF_Score",
                                            "myCAF_Score"))


pdf("./sFeaturePlot_CAF_Score_Stroma.pdf",width = 8,height = 8)

SpatialFeaturePlot(subset_new_sct.3,features = colnames(subset_new_sct.3@meta.data)[c(12:13)],keep.scale = "all",images = c("inmage","inmage.1"),crop = F,pt.size.factor = 7000)

dev.off()

Idents(subset_new_sct.3) <- "orig.ident"

pdf("./VlnPlot_CAF_Score_Stroma.pdf",width = 4,height = 3)
for (i in c(colnames(subset_new_sct.3@meta.data)[c(12:13)])) {
  p1 = VlnPlot(subset_new_sct.3,features = i,idents = c("C15","C16"),cols = c("C15"= "#F8766D" ,"C16" = "#00BFC4"),
               y.max = max(subset_new_sct.3@meta.data[,i])+0.2)+
    stat_compare_means(comparisons = list(c("C15","C16")), label = "p.signif")
  print(p1)
}


dev.off()



subset_new_sct.3 <- subset(new_sct.3,subset = new.cluster.ids %in% c("Immune") & orig.ident %in% c("C15","C16"))
Idents(subset_new_sct.3) <- "orig.ident"


pdf("./sFeaturePlot_ImmuneMarkers_Immune.pdf",width = 6,height = 3)
for (i in c("CD40","CD4","CD8A","LYZ","CD14","CD19")) {
  p1 = SpatialFeaturePlot(subset_new_sct.3,features = i,images = c("inmage","inmage.1"),keep.scale = "all",crop = F,pt.size.factor = 7000)
  print(p1)
}
dev.off()

pdf("./VlnPlot_ImmuneMarkers_Immune.pdf",width = 3,height = 3)

for (i in c("CD40","CD19","CD4","CD8A","CD14","LYZ")) {
  
  p2 = VlnPlot(subset_new_sct.3,features = i,idents = c("C15","C16"),cols = c("C15"= "#F8766D" ,"C16" = "#00BFC4"),
               y.max = max(subset_new_sct.3@assays$Spatial@data[i,])+0.2)+
    stat_compare_means(comparisons = list(c("C15","C16")), label = "p.signif")
  print(p2)
}

dev.off()

```
# 6. spacet draw edge


## 1. define tumor%, test 0.6,0.7,0.8
```{r}
# Identify the Tumor-Stroma Interface 
C15_spacet <- SpaCET.identify.interface(C15_spacet,MalignantCutoff = 0.8)
C16_spacet <- SpaCET.identify.interface(C16_spacet,MalignantCutoff = 0.8)
```

## 2. merge the cell_type percentage data
```{r}
## extract tumor
# tumor <- subset(new_sct.3,subset= new.cluster.ids %in% c("IPMN","PDAC"))
tumor <- new_sct.3

#### extract IPMN spots
seurat <- tumor@images[["inmage"]]
seurat <- seurat@coordinates
seurat$position <- paste0(seurat$row,sep= "x",seurat$col)      
seurat$barcode <- rownames(seurat)

ctdata <- as.data.frame(t(C15_spacet@results$deconvolution$propMat))
ctdata$position <- rownames(ctdata)
pdata <- as.data.frame(cbind(t(C15_spacet@results$CCI$interface),
                             t(C15_spacet@results$CCI$interaction$groupMat)))
pdata$position <- rownames(pdata) 

#merge
tmp_m.pdata <- merge(merge(pdata,seurat,by = "position"),ctdata,by = "position")
#### extract PDAC spots
seurat <- tumor@images[["inmage.1"]]
seurat <- seurat@coordinates
seurat$position <- paste0(seurat$row,sep= "x",seurat$col)      
seurat$barcode <- rownames(seurat)

ctdata <- as.data.frame(t(C16_spacet@results$deconvolution$propMat))
ctdata$position <- rownames(ctdata)
pdata <- as.data.frame(cbind(t(C16_spacet@results$CCI$interface),
                             t(C16_spacet@results$CCI$interaction$groupMat)))
pdata$position <- rownames(pdata) 

#merge
tmp_m.pdata_2 <- merge(merge(pdata,seurat,by = "position"),ctdata,by = "position")


# merge ipmn and pdac

tmp_m.pdata_3 <- rbind(tmp_m.pdata,tmp_m.pdata_2[,c(1:46)])


tmp <- data.frame(barcode = rownames(tumor@meta.data), tumor@meta.data)
tmp <- merge(tmp,tmp_m.pdata_3,by = "barcode")
table(tmp$orig.ident)
## ipmn+pdac=903 spots, only 901 spots have SpaCET data
tumor2 <- tumor[,tmp$barcode]
identical(tmp$barcode,rownames(tumor2@meta.data))
rownames(tmp) <- tmp$barcode
colnames(tmp)
colnames(tumor2@meta.data)
tumor2@meta.data <-tmp
table(tumor2$Interface,tumor2$new.cluster.ids)

SpatialDimPlot(tumor2,group.by = "new.cluster.ids",pt.size.factor = 7000,crop = F,cols = new.cluster.ids_color)
rm(seurat16,tmp,tmp_m.pdata,C16_ctdata,C16_pdata,C15,C16,ctdata,p1,p2,pdata,seurat,tmp_m.pdata_2,tmp_m.pdata_3)

# spacet_06_tumor <- tumor2
spacet_08_tumor <- tumor2

save(spacet_08_tumor,file = "Rdata/spacet_08_tumor.RData")
pdf("SpaCET_pdf/interface.pdf",height = 5,width = 9)
SpatialDimPlot(tumor2,group.by = "Interface",pt.size.factor = 7000,crop = F,cols = c("Interface"= "black",
                                                                                     "Tumor" = "yellow",
                                                                                     "Stroma" = "lightgreen"),alpha = 0.5)
dev.off()


## 3. compare the non-mal cell type within the tumor region
```{r}
data <- spacet_08_tumor@meta.data

sub_data <- subset(data, subset = Interface %in% c("Tumor","Interface"))
sub_data <- sub_data[,c(1,2,5,12,14,21:30,53:55,32:33)]

# Subset data for PDAC and IPMN clusters
pdac_data <- subset(sub_data, orig.ident == "C16")
ipmn_data <- subset(sub_data, orig.ident == "C15")

# Calculate the mean and SD of Malignant, CAF, and Endothelial for each cluster
pdac_means <- colMeans(pdac_data[, c(7:20)])
ipmn_means <- colMeans(ipmn_data[, c(7:20)])

prop <- data.frame(row.names = colnames(pdac_data)[7:20],HG_IPMN = ipmn_means,PDAC = pdac_means)



library(ggplot2)
library(reshape2)
library(scales)
prop2 <- data.frame(celltype = rownames(prop),prop)
prop3 <- melt(prop2)

colors12 <- c("#FF5733", "#8B008B", "#00FFFF", "#FFD700", "#00FF7F",
              "#DC143C", "#00BFFF", "#FF1493", "#1E90FF", "#FF69B4",
              "#FF8C00", "#9400D3")


colors12 <-  c("#F8766D" ,"#00BFC4")

pdf("SpaCET_pdf/within_tumor_08.pdf",width = 8,height = 6)
ggplot(prop3,aes(x = factor(celltype,levels = unique(celltype)),
                 y = ifelse(variable == "PDAC",value,-value),
                 fill = variable))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = colors12)+
  geom_text(aes(label = percent(value,0.01),
                y = ifelse(variable == "PDAC",value,-value),
                hjust = ifelse(variable =="PDAC",0,1)),  ## 标注相对于bar的位置 0和1 刚好在bar上
            size = 4,color = "black")+
  coord_flip(ylim = c(-0.2,0.25))+
  theme_classic() #### eps 700*650
dev.off()

for (i  in 7:20) {
  print(colnames(sub_data)[i])
  print(wilcox.test(sub_data[,i] ~ sub_data$orig.ident, data = sub_data))
}

```

### a. back calculate cell number
```{r}
spot_cell_number = 20
# Subset data for PDAC and IPMN clusters
pdac_data <- subset(sub_data, orig.ident == "C16" & Malignant < 1)
ipmn_data <- subset(sub_data, orig.ident == "C15" & Malignant < 1)

pdac_data[,c(7:20)] <- pdac_data[,c(7:20)] / rowSums(pdac_data[,c(7:20)])
# head(rowSums(pdac_data[,c(7:20)]))
ipmn_data[,c(7:20)] <- ipmn_data[,c(7:20)] / rowSums(ipmn_data[,c(7:20)])

pdac_means <- colMeans(pdac_data[, c(7:20)])
ipmn_means <- colMeans(ipmn_data[, c(7:20)])

prop <- data.frame(row.names = colnames(pdac_data)[7:20],HG_IPMN = ipmn_means,PDAC = pdac_means)



library(ggplot2)
library(reshape2)
library(scales)
prop2 <- data.frame(celltype = rownames(prop),prop)
prop3 <- melt(prop2)

colors12 <-  c("#F8766D" ,"#00BFC4")

pdf("SpaCET_pdf/within_tumor_08.pdf",width = 8,height = 6)
ggplot(prop3,aes(x = factor(celltype,levels = unique(celltype)),
                 y = ifelse(variable == "PDAC",value,-value),
                 fill = variable))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = colors12)+
  geom_text(aes(label = percent(value,0.01),
                y = ifelse(variable == "PDAC",value,-value),
                hjust = ifelse(variable =="PDAC",0,1)),  ## 
            size = 4,color = "black")+
  coord_flip(ylim = c(-0.5,0.5))+
  theme_classic() #### eps 700*650
dev.off()


# 统计占比差异
sub_data2 <- rbind(pdac_data,ipmn_data)
for (i  in 7:20) {
  print(colnames(sub_data2)[i])
  print(wilcox.test(sub_data2[,i] ~ sub_data2$orig.ident, data = sub_data2))
}

```

### b.rm malignant = 1 and plot cell type
```{r}
data <- spacet_08_tumor

sub_data <- subset(data, subset = Interface %in% c("Tumor","Interface"))
sub_data@meta.data <- sub_data@meta.data[,c(1,2,5,12,14,21:30,53:55,32:33)]

sub_data <- subset(sub_data,subset = Malignant != 1)
colnames(sub_data@meta.data)


CAF.gmt <- read.csv("./Genesets/iCAF_myCAF_PMC6727976.csv")

new_sct.3 <- AddModuleScore(new_sct.3,
                            features = list(CAF.gmt[CAF.gmt$term == "iCAF",2],
                                            CAF.gmt[CAF.gmt$term == "myCAF",2]),
                            ctrl  = 100,
                            name = c("iCAF_Score",
                                     "myCAF_Score"))

tmp <- new_sct.3[,rownames(sub_data@meta.data)]
identical(colnames(tmp),rownames(sub_data@meta.data))

sub_data@meta.data <- cbind(sub_data@meta.data,tmp@meta.data[,c(12,13)])


pdf("SpaCET_pdf/within_tumor_08_rm_malignant.pdf",width = 10,height = 7)
for (i in colnames(sub_data@meta.data[7:22]) ) {
  p = SpatialFeaturePlot(sub_data,features = i,crop = F,images = c("inmage","inmage.1"),pt.size.factor = 7000)+NoLegend()
  print(p)
  p1 = VlnPlot(sub_data,features = i,idents = c("C15","C16"),cols = c("C15"= "#F8766D" ,"C16" = "#00BFC4"),
               y.max = max(sub_data@meta.data[,i])+0.2)+
    stat_compare_means(comparisons = list(c("C15","C16")), label = "p.signif")
  print(p1)
}
dev.off()


for (i  in 7:22) {
  print(colnames(sub_data@meta.data)[i])
  print(wilcox.test(sub_data@meta.data[,i] ~ sub_data@meta.data$orig.ident, data = sub_data@meta.data))
}


```



## 4. compare the non-mal cell type outside of the interface
```{r}
sub_data <- subset(data, subset = Interface %in% c("Stroma"))
sub_data <- sub_data[,c(1,2,5,12,14,21:30,53:55,32:33)]

# Subset data for PDAC and IPMN clusters
pdac_data <- subset(sub_data, orig.ident == "C16")
ipmn_data <- subset(sub_data, orig.ident == "C15")

# Calculate the mean and SD of Malignant, CAF, and Endothelial for each cluster
pdac_means <- colMeans(pdac_data[, c(7:20)])
ipmn_means <- colMeans(ipmn_data[, c(7:20)])

prop <- data.frame(row.names = colnames(pdac_data)[7:20],HG_IPMN = ipmn_means,PDAC = pdac_means)

prop2 <- data.frame(celltype = rownames(prop),prop)
prop3 <- melt(prop2)



pdf("SpaCET_pdf/outside_tumor_08.pdf",width = 8,height = 6)
ggplot(prop3,aes(x = factor(celltype,levels = unique(celltype)),
                 y = ifelse(variable == "PDAC",value,-value),
                 fill = variable))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = colors12)+
  geom_text(aes(label = percent(value,0.01),
                y = ifelse(variable == "PDAC",value,-value),
                hjust = ifelse(variable =="PDAC",0,1)),  ## 标注相对于bar的位置 0和1 刚好在bar上
            size = 4,color = "black")+
  coord_flip(ylim = c(-0.2,0.5))+
  theme_classic() 
dev.off()
```



# 7.cnv\monocle for Ductal \> IPMN \> PDAC

```{r}
infercnvApp::infercnvApp()
```

## 1.cnv/monocle3 use xin chuan

```{r}
save(ductal_tumor,file = "ductal_tumor.RData")
```

tried monocle3, but did not have a good result for interpretation

## 2. unclustered

```{r}
infer_res <- infercnv::run(infer,cutoff = 0.1,out_dir = "siCNV/",analysis_mode = "subclusters",cluster_by_groups = F,HMM = F,denoise = F,num_threads = 2)
```

## 3. fine-tune with cluster by group

```{r}
infer_res <- infercnv::run(infer,cutoff = 0.1,out_dir = "siCNV/",analysis_mode = "subclusters",cluster_by_groups = T,HMM = T,denoise = T,num_threads = 2)
```



## 4. manually cluster by subcluster

```{r}
clustering <- read.table("siCNV/infercnv.observation_groupings.txt")

table(clustering$Dendrogram.Group)

```

## 5. replace

```{r}
clustering$Dendrogram.Group[clustering$Dendrogram.Group == "IPMN_s1"] <- "Clone_A"
clustering$Dendrogram.Group[clustering$Dendrogram.Group == "IPMN_s2"] <- "Clone_A"
clustering$Dendrogram.Group[clustering$Dendrogram.Group == "IPMN_s3"] <- "Clone_A"
clustering$Dendrogram.Group[clustering$Dendrogram.Group == "IPMN_s4"] <- "Clone_A"
clustering$Dendrogram.Group[clustering$Dendrogram.Group == "IPMN_s5"] <- "Clone_A"
clustering$Dendrogram.Group[clustering$Dendrogram.Group == "IPMN_s6"] <- "Clone_A"
clustering$Dendrogram.Group[clustering$Dendrogram.Group == "IPMN_s7"] <- "Clone_B"

clustering$Dendrogram.Group[clustering$Dendrogram.Group == "PDAC_s1"] <- "Clone_E"
clustering$Dendrogram.Group[clustering$Dendrogram.Group == "PDAC_s2"] <- "Clone_E"
clustering$Dendrogram.Group[clustering$Dendrogram.Group == "PDAC_s3"] <- "Clone_F"
clustering$Dendrogram.Group[clustering$Dendrogram.Group == "PDAC_s4"] <- "Clone_G"
clustering$Dendrogram.Group[clustering$Dendrogram.Group == "PDAC_s5"] <- "Clone_C"
clustering$Dendrogram.Group[clustering$Dendrogram.Group == "PDAC_s6"] <- "Clone_E"
clustering$Dendrogram.Group[clustering$Dendrogram.Group == "PDAC_s7"] <- "Clone_D"
clustering$Dendrogram.Group[clustering$Dendrogram.Group == "PDAC_s8"] <- "Clone_G"
clustering$Dendrogram.Group[clustering$Dendrogram.Group == "PDAC_s9"] <- "Clone_C" # can not decided by eyes
clustering$Dendrogram.Group[clustering$Dendrogram.Group == "PDAC_s10"] <- "Clone_C" # can not decided by eyes
clustering$Dendrogram.Group[clustering$Dendrogram.Group == "PDAC_s11"] <- "Clone_C" # can not decided by eyes
clustering$Dendrogram.Group[clustering$Dendrogram.Group == "PDAC_s12"] <- "Clone_C" # can not decided by eyes

write.table(clustering,"siCNV/infercnv.observation_groupings.txt",quote = T)

```

## 6. assign clone back to seurat obj

```{r}
ductal <- subset(ductal_tumor,subset = new.cluster.ids == "Ductal")@meta.data

label <- read.table("ductal_tumor_clone.tsv")
colnames(label) <- c("Barcode","Histology")

rownames(label) <- label$Barcode
label <- label[colnames(ductal_tumor),]

identical(rownames(ductal_tumor@meta.data),label$Barcode)

ductal_tumor$clone <- label$Histology

table(ductal_tumor$clone,ductal_tumor$orig.ident)
table(ductal_tumor$clone,ductal_tumor$new.cluster.ids)
table(ductal_tumor$orig.ident,ductal_tumor$new.cluster.ids)

```

### 6.1 heatmap show clone/sample/cluster relationship
```{r}
heatmap_data = ductal_tumor@meta.data[,c(6,11,12)] 
# heatmap_data$new.cluster.ids2 <- ifelse(heatmap_data$label2 == "Invasive_PDAC" & heatmap_data$new.cluster.ids == "IPMN", "PDAC_IPMN",heatmap_data$new.cluster.ids) 

heatmap_data = subset(heatmap_data,clone != "Ductal") 
heatmap_data = heatmap_data %>% arrange(clone,new.cluster.ids)

heatmap_data <- HeatmapAnnotation(PathologyLabel = heatmap_data$label2,
                                  Cluster = heatmap_data$new.cluster.ids, 
                                  Clone = heatmap_data$clone,
                                  col = list(Clone = clone_color,
                                             PathologyLabel = c("HG_IPMN" = "#F8766D" ,"Invasive_PDAC"="#00BFC4"),
                                             Cluster = c("PDAC"="#66C2A5","IPMN"="#A6D854"),
                                             annotation_width = unit(c(1, 4), 'cm'),
                                             gap = unit(1, 'mm')))
pdf("siCNV_pdf/clone_label_cluster.pdf",width = 10,height = 5)
plot(heatmap_data)
dev.off()
```


### 6.2 spatial clone plot
```{r}


SpatialDimPlot(ductal_tumor,group.by = "orig.ident",alpha = 0.5)
DimPlot(ductal_tumor,group.by = "clone")
DimPlot(ductal_tumor,group.by = "new.cluster.ids")
DimPlot(ductal_tumor,group.by = "orig.ident")


# Generate a discrete color palette with 8 colors 
clone_color <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f")

names(clone_color) <- c("Clone_A","Clone_B","Clone_C","Clone_D","Clone_E",
                        "Clone_F","Clone_G","Ductal")

# Plot the colors
par(mar = c(3, 3, 2, 2))
plot(1:7, rep(1, 7), pch = 19, cex = 4, col = clone_color, xlab = "", ylab = "", xaxt = 'n', yaxt = 'n', bty = 'n')+
  text(1:8, rep(1, 8) - 0.05, labels = clone_color, pos = 1, cex = 0.8, col = "black")+
  text(1:8, rep(1, 8) - 0.15, labels = names(clone_color), pos = 1, cex = 0.8, col = "black")


pdf("siCNV_pdf/SDimplot_IPMN_Clone.pdf",width = 8,height = 5)
SpatialDimPlot(ductal_tumor,group.by = "clone",images = "inmage",alpha = 0.8,pt.size.factor = 7000,crop = F,cols = clone_color)
dev.off()
pdf("siCNV_pdf/SDimplot_PDAC_Clone.pdf",width = 8,height = 5)
SpatialDimPlot(ductal_tumor,group.by = "clone",images = "inmage.1",alpha = 0.8,pt.size.factor = 7000,crop = F,cols =clone_color)
dev.off()
write.csv(ductal_tumor@meta.data,file = "siCNV_pdf/clone_distance.csv")


# we do not need to assign clones back to the whole slide

```
### 6.3 clone diff
```{r}
Idents(ductal_tumor) <- "clone"
clone_deg <- FindAllMarkers(subset(ductal_tumor,clone != "Ductal"),only.pos = T)
clone_deg <- clone_deg[clone_deg$p_val_adj<0.05,]
table(clone_deg$cluster)

CGF <- subset(ductal_tumor,clone %in% c("Clone_C","Clone_G","Clone_F") )
CGF <- CreateSeuratObject(counts = CGF@assays$Spatial@counts,meta.data = CGF@meta.data)
CGF <- NormalizeData(CGF)
CGF <- ScaleData(CGF)
Idents(CGF) <- "clone"


DoHeatmap(CGF,features = subset(clone_deg, cluster %in% c("Clone_C","Clone_G","Clone_F"))$gene)

library(enrichR)

dbs <- listEnrichrDbs()
dbs <- c("WikiPathway_2023_Human","MSigDB_Oncogenic_Signatures","Metabolomics_Workbench_Metabolites_2022","MSigDB_Hallmark_2020","KEGG_2021_Human","GO_Molecular_Function_2023","GO_Cellular_Component_2023","GO_Biological_Process_2023","Enrichr_Submissions_TF-Gene_Coocurrence")

enriched <- enrichr(subset(clone_deg, cluster %in% c("Clone_F"))$gene, dbs)

pdf("siCNV_pdf/Barplot_enrichr_Clone_F.pdf",width = 8,height = 8)
for (i in c(1:length(dbs))) {
  
  p1 = plotEnrich(enriched[[i]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",title = paste0(dbs[i]))
  print(p1)
}
dev.off()


```



might do not have so much difference in transcriptomic

## 7.create new infercnv object with clone

```{r}
infer_clone <- infercnv::CreateInfercnvObject(raw_counts_matrix = ensemble,
                                              gene_order_file = "./siCNV_GeneOrderFile.tsv",
                                              annotations_file = "./ductal_tumor_clone.tsv",
                                              delim = "\t",
                                              ref_group_names = "Ductal",
                                              chr_exclude = "chrM")
# perform infercnv operations to reveal cnv signal
infer_res2 <- infercnv::run(infer_clone,cutoff = 0.1,
                            out_dir = "siCNV_clone/",
                            cluster_by_groups = T,
                            HMM = T,
                            denoise = T,
                            num_threads = 2)
```

## 7.1. infercnv plot for pdf
```{r}

plot_cnv(infer_res2,out_dir = "siCNV_pdf/",output_format = "png",output_filename = "infercnv_by_clone_2",png_res = 300,cluster_by_groups = T,contig_cex = 3)


```


## 8.explore genes in cnv heatmap

```{r}
# R libraries required include:
# infercnvNGCHM
# NGCHM
# tsvio
# inferCNV
# RcolorBrewer

library(tsvio)
library(NGCHM)
library(infercnvNGCHM)
infer_res_symbol <- infer_res2
# a <- infer_res_symbol@count.data
a <- infer_res_symbol@expr.data
gene.df <- clusterProfiler::bitr(rownames(a),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db)
gene.df <- gene.df[!duplicated(gene.df$ENSEMBL),]
gene.df <- gene.df[!duplicated(gene.df$SYMBOL),]
#gene.df <- read.table("GeneID.txt",sep = "\t",header = T)
a <- as.matrix(a[gene.df$ENSEMBL,])
class(a)
identical(rownames(a),gene.df$ENSEMBL)


rownames(a)<- gene.df$SYMBOL
a <- as.matrix(as.data.frame(a))

# count.data<- a
expr.data<- a

# infer_res_symbol@count.data <- count.data
infer_res_symbol@expr.data <- expr.data

c <- infer_res_symbol@gene_order
c <- c[rownames(a),]
c$symbol <- gene.df$SYMBOL
identical(rownames(c),gene.df$ENSEMBL)
rownames(c) <- gene.df$SYMBOL

infer_res_symbol@gene_order <- infer_res_symbol@gene_order[,1:3]

# length(rownames(a))
# table(duplicated(gene.df$Gene.stable.ID))
# table(duplicated(gene.df$HGNC.symbol))
# x = gene.df[gene.df$ENSEMBL ==  "ENSG00000183426",]

rm(a,c,count.data,expr.data)
ngchm(infercnv_obj = infer_res_symbol,
      out_dir              = ".",
      path_to_shaidyMapGen = "./ShaidyMapGen.jar",
      gene_symbol          = "bio.gene.hugo")

```

## 9.how to find clone status and clone tree
i6 HMM: a six-state CNV model that predicts the following CNV levels:
  state 1 : 0x = complete loss
state 2 : 0.5x = loss of one copy
state 3 : 1x = neutral
state 4 : 1.5x = addition of one copy
state 5 : 2x = addition of two copies
state 6 : 3x = essentially a placeholder for >2x copies but modeled as 3x.
```{r}
cnv <- read.table("./siCNV_clone/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat",header = T) 

cnv$group <- substr(cnv$cell_group_name,start = 1,stop = 7)

table(cnv$group,cnv$state)

chr_pq <- read.csv("chr_pq_position.csv")

cnv <- subset(cnv, subset = state != 3)
cnv <- subset(cnv,subset = group %in% c("Clone_A","Clone_B","Clone_C","Clone_D","Clone_E","Clone_F","Clone_G"))

cnv$state2 <- cnv$state-3


total_cnv <- cnv %>% group_by(group) %>% summarise(state2_sum = sum(state2, na.rm = TRUE))
write.csv(cnv,file = "siCNV_pdf/cnv.csv")
write.csv(total_cnv,file = "siCNV_pdf/total_cnv.csv")

get_cnv_gain_loss <- function(cnv){
  for (i in 1:nrow(cnv)){
    print(i)
    cnv1 <- cnv[i,]
    chr_pq1 <- chr_pq[chr_pq$chr==cnv1$chr,]
    if(cnv1$end < chr_pq1[1,"end"]) {
      cnv$cnv_type[i] <- paste0(cnv1$chr,"p")
    }else if(cnv1$start > chr_pq1[2,"start"]){
      cnv$cnv_type[i] <- paste0(cnv1$chr,"q")
    }else{
      cnv$cnv_type[i] <- paste0(cnv1$chr,"p,q")
    }
    
    if(cnv1$state > 3){
      cnv$cnv_type[i] <- paste0(cnv$cnv_type[i],"_gain")
    }else if(cnv1$state < 3){
      cnv$cnv_type[i] <- paste0(cnv$cnv_type[i],"_loss")
    }
    
    # if(grepl(",",cnv$cnv_type[i])){
    #   cha1 <- gsub("p,q.*","",cnv$cnv_type[i])
    #   cha2 <- gsub(".*_","",cnv$cnv_type[i])
    #   cnv$cnv_type[i] <- paste0(cha1,"q_",cha2,",",cha1,"p_",cha2)
  }
  
  cnv1 <- unique(cnv[,c(7,8)])
  # cnv1.1 <- cnv1[!grepl(",",cnv1$cnv_type),]
  # cnv1.2 <- cnv1[grepl(",",cnv1$cnv_type),]
  # cnv1.21 <- cnv1.2
  # cnv1.21$cnv_type <- gsub(",.*","",cnv1.21$cnv_type)
  # cnv1.22 <- cnv1.2
  # cnv1.22$cnv_type <- gsub(".*,","",cnv1.22$cnv_type)
  # cnv1 <- unique(do.call(rbind,list(cnv1.1,cnv1.21,cnv1.22)))
  cnv1 <- unique(do.call(rbind,list(cnv1)))
  table(cnv1$group)
  cnv2 <- reshape2::dcast(cnv1,cnv_type~group)
  index <- as.numeric(gsub("._.*|chr","",cnv2$cnv_type))
  cnv2 <- cnv2[order(index),]
  return(cnv2)
}


cnv2 <- get_cnv_gain_loss(cnv)

write.csv(cnv2,file = "./cnv_gain_loss.csv")
```

from the cnv plot
Clone_A = IPMN
Clone_B/D/G shared with chr5 amplification
Clone_F has exclusive chr2 amplification

compare BDG/A with chr5 high
compare F/A with chr2 high

```{r}
Idents(ductal_tumor) <- "clone"
deg_BDG <- FindMarkers(ductal_tumor,ident.1 = c("Clone_B","Clone_D","Clone_G"),ident.2 = "Clone_A")

gene.chr <- read.table("siCNV_clone/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat",header = T)

deg_BDG$symbol <- rownames(deg_BDG)
ensembl <- bitr(deg_BDG$symbol,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = org.Hs.eg.db,drop = T)
deg_BDG <- deg_BDG[ensembl$SYMBOL,]
deg_BDG$ensembl <- ensembl$ENSEMBL
deg_BDG <- merge(deg_BDG,by.x = "ensembl",gene.chr, by.y = "gene")

deg_BDG$group <- substr(deg_BDG$cell_group_name ,1,7)

deg_BDG <- clusterProfiler::filter(deg_BDG,group %in% c("Clone_B","Clone_D","Clone_G") & avg_log2FC > 0 & p_val_adj < 0.05 & chr== "chr5")

write.csv(deg_BDG,file = "deg_BDGvsA_chr5.csv")


deg_F <- FindMarkers(ductal_tumor,ident.1 = c("Clone_F"),ident.2 = "Clone_A")
deg_F$symbol <- rownames(deg_F)
ensembl <- bitr(deg_F$symbol,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = org.Hs.eg.db,drop = T)
deg_F <- deg_F[ensembl$SYMBOL,]
deg_F$ensembl <- ensembl$ENSEMBL
deg_F <- merge(deg_F,by.x = "ensembl",gene.chr, by.y = "gene")

deg_F$group <- substr(deg_F$cell_group_name ,1,7)

deg_F <- clusterProfiler::filter(deg_F,group %in% c("Clone_F") & avg_log2FC > 0 & p_val_adj < 0.05 & chr =="chr2")
table(deg_F$chr)
write.csv(deg_F,file = "deg_FvsA_chr2.csv")

```


# 8. CellChat

```{r}
library(patchwork) 
library(ggplot2)
library(CellChat)
library(ggalluvial)
library(svglite)
```

## 1. subset for each sample (remove the normal spots)
```{r}
a = subset(new_sct.3,subset = orig.ident == "C15")
a = subset(new_sct.3,subset = new.cluster.ids %in% c("IPMN","Stroma","Immune"))
cellchat_15 <- createCellChat(object = a@assays$Spatial@data,meta = a@meta.data,group.by = "new.cluster.ids")
a = subset(new_sct.3,subset = orig.ident == "C16")
a = subset(new_sct.3,subset = new.cluster.ids %in% c("PDAC","IPMN","Stroma","Immune"))
cellchat_16 <- createCellChat(object = a@assays$Spatial@data,meta = a@meta.data,group.by = "new.cluster.ids")
rm(a)
showDatabaseCategory(CellChatDB.human)
dplyr::glimpse(CellChatDB.human$interaction)
```
## 2. process the cellchat using just secreated signaling
```{r}
cellchat_obj <- cellchat_15

# cellchat_obj <- cellchat_16
cellchat_obj@idents = droplevels(cellchat_obj@idents, exclude = setdiff(levels(cellchat_obj@idents),unique(cellchat_obj@idents)))


## use all database
cellchat_obj@DB <- subsetDB(CellChatDB.human,search = "Secreted Signaling")
## pre process
cellchat_obj <- subsetData(cellchat_obj)
cellchat_obj <- cellchat_obj %>% identifyOverExpressedGenes() %>% identifyOverExpressedInteractions() %>% smoothData(adj = PPI.human)
## infer cellchat network
cellchat_obj <- computeCommunProb(cellchat_obj,raw.use = FALSE)
cellchat_obj <- filterCommunication(cellchat_obj,min.cells = 10)
cellchat_obj <- cellchat_obj %>% computeCommunProbPathway() %>% aggregateNet()

cellchat_15 <- cellchat_obj
# cellchat_16 <- cellchat_obj

save(cellchat_15,cellchat_16,file = "./Rdata/cellchat.Rdata")

```


## 3. find best K (skip)
```{r}
## analysis network find KKKKK
cellchat_obj <- netAnalysis_computeCentrality(cellchat_obj)
selectK(cellchat_obj,pattern = "outgoing")
cellchat_obj <- identifyCommunicationPatterns(cellchat_obj,pattern = "outgoing",k = 3)


##
df.net <- subsetCommunication(cellchat_obj)
table(df.net$pathway_name)
hist(df.net$prob)
df.net_filter <- filter(df.net,df.net$source %in% c("Immune","Stroma") &
                          df.net$target == "PDAC" & df.net$prob > 0.02)


groupSize <- as.numeric(table(cellchat_obj@idents))
netVisual_circle(cellchat_obj@net$count,vertex.weight = groupSize,weight.scale = T,label.edge = F,title.name = "Number of Interactions")
netVisual_circle(cellchat_obj@net$weight,vertex.weight = groupSize,weight.scale = T,label.edge = F,title.name = "Weight of Interactions")


## contribution
pathway = "SPP1"
netAnalysis_contribution(cellchat_obj, signaling = pathway)



vertex.receiver = seq(1,4)
netVisual_aggregate(cellchat_obj, signaling = pathway,  vertex.receiver = vertex.receiver)

par(mfrow=c(1,1))
netVisual_aggregate(cellchat_obj, signaling = pathway, layout = "circle")


# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_obj, signaling = pathway, color.heatmap = "Reds")

netVisual_bubble(cellchat_obj,sources.use = "Immune",targets.use = "PDAC",remove.isolate = F)
netVisual_chord_gene(cellchat_obj,sources.use = "Immune",targets.use = "PDAC",lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat_obj,sources.use = "Immune",targets.use = "PDAC",lab.cex = 0.5,slot.name = "netP")


```



## 4. merge
```{r}

group.new = levels(cellchat_16@idents)
cellchat_15<- liftCellChat(cellchat_15,group.new = group.new)

object <- list(IPMN = cellchat_15,PDAC = cellchat_16)
merge_cellchat <- mergeCellChat(object , add.names = names(object),cell.prefix = T)

pdf("./cellchat/BarPlot_primary_comparison.pdf",width = 6,height = 4)
compareInteractions(merge_cellchat, show.legend = F, group = c(1,2)) | compareInteractions(merge_cellchat, show.legend = F, group = c(1,2), measure = "weight") 
dev.off()
netAnalysis_signalingChanges_scatter(merge_cellchat, idents.use = "Stroma")

## change the IPMN and PDAC name to Tumor
merge_cellchat <- computeNetSimilarityPairwise(merge_cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
merge_cellchat <- netEmbedding(merge_cellchat, type = "functional",umap.method = "uwot")
#> Manifold learning of the signaling networks for datasets 1 2
merge_cellchat <- netClustering(merge_cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2

# Visualization in 2D-space
pdf("./cellchat/Clusters.pdf",width = 10,height = 6)
netVisual_embeddingPairwise(merge_cellchat, type = "functional", label.size = 3.5)
dev.off()
#> 2D visualization of signaling networks from datasets 1 2

df.net <- subsetCommunication(merge_cellchat)

pdf("./cellchat/BarPlot_stacked_info_flow_Stroma.pdf",width = 6,height = 6)
rankNet(merge_cellchat,mode = "comparison",sources.use = c(2),targets.use = c(1, 4) ,stacked = T,do.flip =T, do.stat = T,cutoff.pvalue = 0.05)
dev.off()

pdf("./cellchat/BarPlot_stacked_info_flow_Immune.pdf",width = 6,height = 6)
rankNet(merge_cellchat,mode = "comparison",sources.use = c(3),targets.use = c(1, 4) ,stacked = T,do.flip =T, do.stat = T,cutoff.pvalue = 0.05)
dev.off()

pdf("./cellchat/DotPlot_0.05_Stroma_PDAC_IPMN.pdf",width = 5,height = 12)
netVisual_bubble(merge_cellchat,comparison = c(1, 2),sources.use = c(2),targets.use = c(1, 4) , angle.x = 45,thresh = 0.05,remove.isolate = T,signaling = c("SPP1","PERIOSTIN","PLAU","ANGPTL","TWEAK","BMP","HGF"))
dev.off()
pdf("./cellchat/DotPlot_Immune_PDAC_IPMN.pdf",width = 5,height = 12)
netVisual_bubble(merge_cellchat,comparison = c(1, 2),sources.use = c(3),targets.use = c(1, 4) , angle.x = 45,thresh = 0.05,remove.isolate = T,signaling = c())
dev.off()

pdf("./cellchat/BarPlot_0.05_Stacked_StromaImmune_PDAC_IPMN.pdf",width = 5,height = 6)
rankNet(merge_cellchat,mode = "comparison",sources.use = c(2,3),targets.use = c(1,4) ,stacked = T,do.flip =T, do.stat = T,cutoff.pvalue = 0.05)
dev.off()
pdf("./cellchat/DotPlot_0.05_StromaImmune_PDAC_IPMN.pdf",width = 5,height = 5)
netVisual_bubble(merge_cellchat,comparison = c(1, 2),sources.use = c(2,3),targets.use = c(1, 4) , angle.x = 45,thresh = 0.05,remove.isolate = T,signaling = c("SPP1","PERIOSTIN","PLAU","BMP","HGF","TGFb"))
dev.off()

```





