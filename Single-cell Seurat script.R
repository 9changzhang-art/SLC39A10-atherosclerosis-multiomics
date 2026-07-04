
library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
library(plyr)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(SCpubr)
library(UCell)
library(irGSEA)
library(GSVA)
library(GSEABase)
library(harmony)
library(randomcoloR)
library(CellChat)
library(ggpubr)

setwd("E:\\xin21\\SLC\\dan")
af=readRDS("af8.rds")
afgs="zinc.txt"


pdf(file = "01.vlnplot.pdf",width = 8,height = 5)
VlnPlot(af, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 6)+scale_fill_manual(values = c("#C77CFF","#7CAE00","#00BFC4","#F8766D","#FFD700","#00FF00"))
dev.off()


plot1 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "percent.mt")+ RotatedAxis()
plot2 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "percent.rb")+ RotatedAxis()
plot3 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ RotatedAxis()

pdf(file = "01.corqc.pdf",width =12,height = 5)
plot1+plot2+plot3+plot_layout(ncol = 3)      #plot_layout，patchwork函数，指定一行有几个图片
dev.off()


af <- NormalizeData(af, normalization.method = "LogNormalize", scale.factor = 10000)


af <- FindVariableFeatures(af, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(af), 10)
top10


plot1 <- VariableFeaturePlot(af)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf(file = "01.topgene.pdf",width =7,height = 6)
plot2                   
dev.off()


all.genes <- rownames(af)
af <- ScaleData(af, features = all.genes)


af <- Seurat::RunPCA(af, features = VariableFeatures(object = af))
af <- Seurat::RunTSNE(af,dims = 1:20)
pdf(file = "02.rawtsne.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "tsne",pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right") 
dev.off()
pdf(file = "02.rawpca.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "pca",pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
colaa=distinctColorPalette(100)
pdf(file = "02.raw.tsne.split.pdf",width =8,height =5)
do_DimPlot(sample = af,
           plot.title = "",
           reduction = "tsne",
           legend.position = "bottom",
           dims = c(1,2),split.by = "Type",pt.size =0.5
) 
dev.off()


af <- RunHarmony(af, group.by.vars = "Type")


pdf(file = "03.harmony.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "harmony",pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
af <- Seurat::RunTSNE(af,dims = 1:20,reduction ='harmony')
pdf(file = "03.tsne.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "tsne",pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
collist=c(ggsci::pal_nejm()(8))
names(collist)=names(table(af$Type))
pdf(file = "03.tsne.split.pdf",width =12,height = 7.5)
do_DimPlot(sample = af,
           plot.title = "",
           reduction = "tsne",
           legend.position = "bottom",
           dims = c(1,2),split.by = "Type",pt.size =0.5
) 
dev.off()


collist=c(ggsci::pal_nejm()(8))
names(collist)=names(table(af$Type))

VizDimLoadings(af, dims = 1:2, reduction = "pca")

pdf(file = "04.pc_heatmap.pdf",width =7.5,height = 9)
DimHeatmap(af, dims = 1:20, cells = 1000, balanced = TRUE)
dev.off()

af <- JackStraw(af, num.replicate = 100)
af <- ScoreJackStraw(af, dims = 1:20)
pdf(file = "04.jackstrawplot.pdf",width =7.5,height = 5.5)
JackStrawPlot(af, dims = 1:20)
dev.off()
pdf(file = "04.ElbowPlot.pdf",width =5,height = 4)
ElbowPlot(af,ndims = 30,reduction = "harmony")
dev.off()


afPC=15


af=FindNeighbors(af, dims = 1:afPC, reduction = "harmony")
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1,1.2,1.5,2,2.5,3)) {
  af=FindClusters(af, graph.name = "RNA_snn", resolution = res, algorithm = 1)}
apply(af@meta.data[,grep("RNA_snn_res",colnames(af@meta.data))],2,table)

p2_tree=clustree(af@meta.data, prefix = "RNA_snn_res.")
pdf(file = "04.clustertree.pdf",width =12,height =10)
p2_tree
dev.off()

af=FindNeighbors(af, dims = 1:afPC, reduction = "harmony")
af <- FindClusters(af, resolution = 0.3) #关键！


head(Idents(af), 5)


head(af@meta.data)
table(af@meta.data$seurat_clusters)

af.markers <- FindAllMarkers(af, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(af.markers,file = "05.cluster_markers.csv")


af <- RunUMAP(af, dims = 1:afPC, reduction = "harmony")
af <- RunTSNE(af, dims = 1:afPC, reduction = "harmony")


pdf(file = "05-cluster.UMAP.pdf",width =7,height = 5.5)
DimPlot(af, reduction = "umap", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "05-cluster.TSEN.pdf",width =7,height = 5.5)
DimPlot(af, reduction = "tsne", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()



afgenes=read.table("ppi.hub.txt",header = F,sep = "\t")[,1]
afgss=as.character(read.table(afgs,header = F,sep = "\t")[,1])
af <- AddModuleScore(
  object = af,
  features =list(afgss[afgss %in% rownames(af)]) ,
  ctrl = 100, 
  name = gsub(".txt","",afgs)
)
colnames(af@meta.data)[length(colnames(af@meta.data))] <- gsub(".txt","",afgs) 
colsa = distinctColorPalette(100) 


source(file = "vnplot.R")
gene_sig <- gsub(".txt","",afgs)
comparisons <- list(names(table(af$Type)))
afvp(af,gene_signature = gene_sig, file_name = "06-group.GS_VlnPlot", test_sign = comparisons,pta=0.1,cols=colsa,label="p.format",group="Type",widplot=6,heiplot=5.5) #

comparisons <- list(names(table(af$Type)))
afvp(af,gene_signature = intersect(afgenes,rownames(af)), file_name = "06-group.hub_VlnPlot", test_sign = comparisons,pta=0.1,cols=colsa,label="p.signif",group="Type",widplot=15,heiplot=12,ak=0.9) 

comparisons <- list()
comp=combn(names(table(af$seurat_clusters)),2)
names(table(af$seurat_clusters))
for(j in 1:ncol(comp)){comparisons[[j]]<-comp[,j]}
afvp(af,gene_signature = gene_sig, file_name = "06-cluster.GS_VlnPlot", test_sign = comparisons,pta=0.1,cols=colsa,label="p.signif",group="seurat_clusters",widplot=20,heiplot=8,ak=0.01,split = "Type")

pdf(file = "06-cluster.hub_VlnPlot.pdf",width =10,height = 6)
VlnPlot(af, features = afgenes,group.by = "seurat_clusters", stack=TRUE,cols = colsa, slot = "data")+ NoLegend()   
dev.off()




af <- FindClusters(af, resolution = 0.3)
genes <- list("Smooth muscle cells" = c("ACTA2","MYH11","CNN1"),
              "Macrophages" = c("CD163", "C1QA","C1QB"),
              "T cells" = c("CD3D", "CD3E","TRAC"),
              "B cells" = c("CD79A", "MS4A1", "CD37"),
              "Mast cells" = c("CPA3", "MS4A2", "TPSAB1"),
              "Progenitor_cells" = c("MKI67","ASPM","TOP2A"),
              "Endothelial cells" = c("PECAM1", "VWF", "EGFL7"),
              "Plasma cells"=c("MZB1","JCHAIN","IGLC2"),
              "NK cells" = c("KLRC1", "TRDC", "XCL2"),
              "Fibroblasts" = c("FBLN1", "SFRP2","DCN"),
              "Dendritic cells" = c("CD1C","LGALS2","CLEC10A"),
              "Neutrophils" = c("S100A9", "S100A8","TREM1"),
              "SEM" = c("VCAM1", "COL1A1", "COL1A2", "LGALS3", "FN1")
)
pdf(file = "0850001.ann_cluster_marker.pdf",width =20,height = 7)
do_DotPlot(sample = af,features = genes,dot.scale = 12,sequential.palette = "YlOrRd",legend.length = 50,
           legend.framewidth = 2, font.size =12)
dev.off()



table(af@active.ident)
ann.ids <- c("T cells",  #cluster0
             "Smooth muscle cells",  #cluster1
             "T cells",     
             "Endothelial cells",
             "Macrophages",
             "Neutrophils",
             "Fibroblasts",
             "NK cells",
             "B cells",
             "Progenitor cells",
             "Dendritic cells",#cluster10
             "Mast cells",
             "Plasma cells",
             "SEM",
             "Plasma cells",
             "Endothelial cells",
             "Smooth muscle cells"
)
length(ann.ids)

afidens=mapvalues(Idents(af), from = levels(Idents(af)), to = ann.ids)
Idents(af)=afidens
af$cellType=Idents(af)


pdf(file = "08-ann.scRNA.UMAP.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "umap", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "08-ann.scRNA.TSEN.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "tsne", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()


af.markers <- FindAllMarkers(af, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(af.markers,file = "08.cell_markers.csv")

top5af.markers <- af.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)


pdf(file = "09-cell_marker.hetmap.pdf",width =15,height = 10)
DoHeatmap(af,features = top5af.markers$gene,
          group.colors = colsa) +
  ggsci::scale_colour_npg() +
  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',
                       name = 'Z-score')
dev.off()



colaa=distinctColorPalette(100)
afgenes=read.table("ppi.hub.txt",header = F,sep = "\t")[,1]
pdf(file = "09-cell_FeaturePlot.pdf",width =12,height = 10)
FeaturePlot(af, features = afgenes, cols = c("grey", "red"),min.cutoff = 0.1, max.cutoff = 1,ncol=4,pt.size = 0.5, slot = "counts")    
pdf(file = "09-cell_VlnPlot.pdf",width =8,height = 5)
VlnPlot(af, features = afgenes,group.by = "cellType", stack=TRUE,cols = colaa, slot = "counts")+ NoLegend()   
dev.off()
pdf(file = "09-cell_Group_VlnPlot.pdf",width =8,height = 5)
VlnPlot(af, features = afgenes,group.by = "cellType", stack=TRUE,cols = c("red3","blue3"), slot = "counts", split.by ="Type")
dev.off()
colsa = distinctColorPalette(100)
pdf(file = "09-cell_GS_VlnPlot.pdf",width =8,height = 5)
VlnPlot(af, features = gsub(".txt","",afgs),group.by = "cellType",cols=colsa, split.by ="Type")
dev.off()


source(file = "vnplot.R")
gene_sig <- gsub(".txt","",afgs)

comparisons <- list()
comp=combn(names(table(af$cellType)),2)
names(table(af$cellType))
for(j in 1:ncol(comp)){comparisons[[j]]<-comp[,j]}
afvp(af,gene_signature = gene_sig, file_name = "09-cell.GS.stat_VlnPlot", test_sign = comparisons,pta=0.1,cols=colsa,label="p.signif",group="cellType",widplot=10,heiplot=10,ak=0.9)

comparisons <- list()
comp=combn(names(table(af$cellType)),2)
names(table(af$cellType))
for(j in 1:ncol(comp)){comparisons[[j]]<-comp[,j]}
afvp(af,gene_signature = intersect(afgenes,rownames(af)), file_name = "09-cell.hub.stat_VlnPlot", test_sign = comparisons,pta=0.1,cols=colsa,label="p.signif",group="cellType",widplot=20,heiplot=15,ak=5)


geneselect="SLC39A10"
cellselect="SEM"
af$geneType=ifelse(as.numeric(as.matrix(af@assays$RNA@scale.data)[geneselect,])>median(sort(as.numeric(as.matrix(af[,which(af$cellType %in% c(cellselect))]@assays$RNA@scale.data)[geneselect,]))),paste0("High ",geneselect," ",cellselect),paste0("Low ",geneselect," ",cellselect))
Cellratio <- prop.table(table( af[,which(af$cellType %in% c(cellselect))]$geneType,af[,which(af$cellType %in% c(cellselect))]$Type), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio)[1]="Celltype"
colourCount = length(unique(Cellratio$Celltype))
colaa=distinctColorPalette(100)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Celltype),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Type',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")+   
  scale_fill_manual(values=colaa)
ggsave("102-cell_geneselect_ration.pdf",width = 6,height = 3.5)

afvp(af,gene_signature = geneselect, file_name = "06-cell.geneselect._VlnPlot", test_sign = comparisons,pta=0.1,cols=colsa,label="p.signif",group="Type",widplot=12,heiplot=8,ak=0.9,split = "cellType")

afc=af[,which(af$cellType %in% c(cellselect))]
Idents(afc)=afc$geneType

pdf(file = "102-sg.scRNA.UMAP.pdf",width =7.5,height = 5.5)
DimPlot(afc, reduction = "umap", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "102-sg.scRNA.TSEN.pdf",width =7.5,height = 5.5)
DimPlot(afc, reduction = "tsne", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()


af.markers <- FindAllMarkers(afc, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(af.markers,file = "113.geneGroup_Diff.csv")

library(clusterProfiler)
library(org.Hs.eg.db)
ids=bitr(af.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
af.markers=merge(af.markers,ids,by.x='gene',by.y='SYMBOL')
View(af.markers)

gcSample=split(af.markers$ENTREZID, af.markers$cluster) 

xx <- compareCluster(gcSample,
                     fun = "enrichKEGG",
                     organism = "hsa",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05
)
write.csv(as.data.frame(xx),file = "11.geneGroup_KEGG.csv")
p <- dotplot(xx)
pdf(file = "11-geneGroup_KEGG.pdf",width =8,height = 8)
p +scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))+ theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
dev.off()

xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     #ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05
)
write.csv(as.data.frame(xx),file = "11.geneGroup_GO.csv")
p <- dotplot(xx)
pdf(file = "11-geneGroup_GO.pdf",width =8,height = 8)
p+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60)) + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
dev.off()


