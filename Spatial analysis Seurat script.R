library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Rfast2)

sampleName="GC1"       
setwd("C:\\Users\\huawei\\Desktop\\xinde\\新\\kongjian")     

data = Read10X(".")
image2 <- Read10X_Image(image.dir = "spatial", filter.matrix = TRUE)


stData <- CreateSeuratObject(counts=data, assay="Spatial", project=sampleName)
image2 <- image2[Cells(x = stData)]
DefaultAssay(stData = image2) <- "Spatial"
stData[["slice1"]] <- image2


plot1 <- VlnPlot(stData, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(stData, features = "nCount_Spatial") + theme(legend.position = "right")
pdf(file="ST01.FeaturePlot.pdf", width=9, height=6)
wrap_plots(plot1, plot2)
dev.off()


stData <- SCTransform(stData, assay = "Spatial", verbose = FALSE)


stData <- RunPCA(stData, assay = "SCT", verbose = FALSE)
stData <- FindNeighbors(stData, reduction = "pca", dims = 1:5)
stData <- FindClusters(stData, verbose = FALSE)
stData <- RunUMAP(stData, reduction = "pca", dims = 1:5)


p1 <- DimPlot(stData, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(stData, label = TRUE, label.size = 3)
pdf(file="ST02.SpatialDimPlot.pdf", width=12, height=6)
p1 + p2
dev.off()


pdf(file="ST03.singleDimPlot.pdf", width=12, height=8)
SpatialDimPlot(stData, cells.highlight = CellsByIdentities(object = stData,
    idents = levels(stData$seurat_clusters)), facet.highlight = TRUE, ncol = 4)
dev.off()


save(stData, file="ST.Rdata")


logFCfilter=1           
adjPvalFilter=0.05     
de_markers <- FindAllMarkers(object = stData,
                   only.pos = FALSE,
                   min.pct = 0.1,
                   logfc.threshold = logFCfilter)
sig.markers=de_markers[(abs(as.numeric(as.vector(de_markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(de_markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="ST04.clusterMarkers.txt",sep="	",row.names=F,quote=F)



stData <- FindSpatiallyVariableFeatures(stData, assay = "SCT", 
			features = VariableFeatures(stData)[1:100],
			selection.method = "moransi")

top.features=rownames(dplyr::slice_min(stData[["SCT"]]@meta.features,
					  moransi.spatially.variable.rank, n = 10))
pdf(file="ST076.SpatialGene.pdf", width=15, height=10)
SpatialFeaturePlot(stData, features = top.features, ncol=3, alpha=c(0.1, 1))
dev.off()



zinc_genes <- c("SLC39A1", "SLC39A2", "SLC39A3", "SLC39A4", "SLC39A5", 
                "SLC39A6", "SLC39A7", "SLC39A8", "SLC39A9", "SLC39A11",
                "SLC39A12", "SLC39A13", "SLC39A14", "SLC30A1", "SLC30A2",
                "SLC30A3", "SLC30A4", "SLC30A5", "SLC30A6", "SLC30A7",
                "SLC30A8", "SLC30A9", "SLC30A10", "MT1A", "MT1B", "MT1E",
                "MT1F", "MT1G", "MT1H", "MT1M", "MT1X", "MT2A", "MT3", "MT4")


all_genes <- rownames(stData)
zinc_genes_present <- intersect(zinc_genes, all_genes)


sink("ST09.ZincGeneInfo.txt")

sct_data <- GetAssayData(stData, assay = "SCT", layer = "data")

if(length(zinc_genes_present) >= 1) {
  zinc_expr <- sct_data[zinc_genes_present, , drop = FALSE]
  if(nrow(zinc_expr) == 1) {
    stData$Zinc_Score <- as.numeric(zinc_expr[1, ])
  } else {
    stData$Zinc_Score <- colMeans(zinc_expr, na.rm = TRUE)
  }
} else {
  stData$Zinc_Score <- 0
  }

zn_q25 <- quantile(stData$Zinc_Score, probs = 0.25, na.rm = TRUE)
zn_q75 <- quantile(stData$Zinc_Score, probs = 0.75, na.rm = TRUE)

stData$Zn_Level <- ifelse(stData$Zinc_Score >= zn_q75, "Zn-high",
                   ifelse(stData$Zinc_Score <= zn_q25, "Zn-low", "Zn-mid"))


p_zinc <- SpatialFeaturePlot(stData, features = "Zinc_Score", 
                              alpha = c(0.1, 1), pt.size.factor = 1.5) +
  scale_fill_gradientn(colors = c("darkblue", "blue", "lightblue", "white", "yellow", "orange", "red"),
                       name = "Zinc Score") +
  ggtitle("Zinc Ion Distribution") +
  theme(legend.position = "right")

pdf(file="ST09.ZincSpatialDistribution.pdf", width=10, height=8)
print(p_zinc)
dev.off()


p_zn_level <- SpatialDimPlot(stData, group.by = "Zn_Level", label = FALSE, pt.size.factor = 1.5) +
  scale_fill_manual(values = c("Zn-high" = "#E74C3C", 
                               "Zn-low" = "#3498DB", 
                               "Zn-mid" = "#95A5A6"),
                    name = "Zn Level") +
  ggtitle("Zinc Ion Level Regions") +
  theme(legend.position = "right")

pdf(file="ST09.ZincLevelRegions.pdf", width=10, height=8)
print(p_zn_level)
dev.off()


p_umap_zinc <- FeaturePlot(stData, features = "Zinc_Score", reduction = "umap") +
  scale_color_gradientn(colors = c("darkblue", "blue", "lightblue", "white", "yellow", "orange", "red"))

pdf(file="ST09.ZincUMAP.pdf", width=8, height=6)
print(p_umap_zinc)
dev.off()

p_vln_zinc <- VlnPlot(stData, features = "Zinc_Score", pt.size = 0.1) + 
  ggtitle("Zinc Score by Cluster")

pdf(file="ST09.ZincScoreVlnPlot.pdf", width=10, height=6)
print(p_vln_zinc)
dev.off()


n_zn_high <- sum(stData$Zn_Level == "Zn-high")
n_zn_low <- sum(stData$Zn_Level == "Zn-low")

if(n_zn_high >= 3 & n_zn_low >= 3) {
  de_zn <- FindMarkers(stData, ident.1 = "Zn-high", ident.2 = "Zn-low", 
                        group.by = "Zn_Level", assay = "SCT")
  de_zn$gene <- rownames(de_zn)
  sig_de_zn <- de_zn[de_zn$p_val_adj < 0.05 & abs(de_zn$avg_log2FC) > 0.5, ]
  write.table(sig_de_zn, file="ST09.ZnHighVsLow_DEG.txt", sep="	", row.names=F, quote=F)

  if(nrow(sig_de_zn) >= 5) {
    top_genes <- head(sig_de_zn$gene[order(abs(sig_de_zn$avg_log2FC), decreasing = TRUE)], 20)
    p_heatmap <- DoHeatmap(subset(stData, Zn_Level %in% c("Zn-high", "Zn-low")), 
                           features = top_genes, group.by = "Zn_Level") +
      ggtitle("Zn-high vs Zn-low DEGs")
    pdf(file="ST09.ZnDEG_Heatmap.pdf", width=12, height=10)
    print(p_heatmap)
    dev.off()
  }
}


save(stData, file="ST_with_ZincAnalysis.Rdata")




ecm_candidate_genes <- c("MMP2", "MMP9", "MMP12", "MMP13", "MMP14", "MMP16",
                        "CTSK", "CTSB", "CTSD", "COL1A1", "COL3A1", "ELN", "FBN1",
                        "TIMP1", "TIMP2", "TIMP3", "LOX", "LOXL1", "LOXL2")

ecm_present <- intersect(ecm_candidate_genes, rownames(stData))


sink("ST10.ECM_Gene_Check.txt")


if(length(ecm_present) > 0) {
  cat("

")
  for(g in ecm_present) {
    expr_vals <- as.numeric(sct_data[g, ])
    cat(g, ": mean=", round(mean(expr_vals), 4), 
        ", sd=", round(sd(expr_vals), 4), 
        ", max=", round(max(expr_vals), 4), 
        ", non-zero=", sum(expr_vals > 0), "/", length(expr_vals), "
")
  }
}
sink()


if(length(ecm_present) >= 1) {
  expr <- sct_data[ecm_present, , drop = FALSE]
  if(nrow(expr) == 1) {
    stData$ECM_Degradation_Score <- as.numeric(expr[1, ])
  } else {
    stData$ECM_Degradation_Score <- colMeans(expr, na.rm = TRUE)
  }
  cat("ECM评分计算完成，使用", length(ecm_present), "个基因
")
} else {
  stData$ECM_Degradation_Score <- 0
  cat("警告：未找到ECM相关基因，ECM评分设为0
")
}

if(length(ecm_present) > 0) {
  gene_vars <- apply(sct_data[ecm_present, , drop=FALSE], 1, var, na.rm=TRUE)
  top_var_gene <- names(which.max(gene_vars))
  stData$Top_ECM_Gene <- as.numeric(sct_data[top_var_gene, ])
  cat("变异最大的ECM基因:", top_var_gene, "，方差=", round(max(gene_vars), 6), "
")
} else {
  top_var_gene <- NULL
}


if(!is.null(top_var_gene)) {
  cor_zn_ecm_module <- cor.test(stData$Zinc_Score, stData$ECM_Degradation_Score, method = "spearman")
  cor_zn_ecm_top <- cor.test(stData$Zinc_Score, stData$Top_ECM_Gene, method = "spearman")

  sink("ST10.Zn_ECM_Correlation.txt")
  cat("=== 锌-ECM相关性分析 ===

")
  cat("【模块评分相关性】
")
  cat("  Spearman rho:", round(cor_zn_ecm_module$estimate, 4), "
")
  cat("  P-value:", format(cor_zn_ecm_module$p.value, digits=4, scientific=TRUE), "

")
  cat("【Top变异基因", top_var_gene, "相关性】
")
  cat("  Spearman rho:", round(cor_zn_ecm_top$estimate, 4), "
")
  cat("  P-value:", format(cor_zn_ecm_top$p.value, digits=4, scientific=TRUE), "
")
  sink()
}


smc_contractile <- c("ACTA2", "MYH11", "TAGLN", "CNN1", "TPM2")
smc_synthetic <- c("S100A4", "VIM", "FN1", "LGALS3", "CD44", "VCAM1")

smc_con_present <- intersect(smc_contractile, rownames(stData))
smc_syn_present <- intersect(smc_synthetic, rownames(stData))

if(length(smc_con_present) > 0) {
  expr <- sct_data[smc_con_present, , drop = FALSE]
  stData$SMC_Contractile_Score <- ifelse(nrow(expr)==1, as.numeric(expr[1,]), colMeans(expr, na.rm=TRUE))
}
if(length(smc_syn_present) > 0) {
  expr <- sct_data[smc_syn_present, , drop = FALSE]
  stData$SMC_Synthetic_Score <- ifelse(nrow(expr)==1, as.numeric(expr[1,]), colMeans(expr, na.rm=TRUE))
}


stData$is_SMC <- FALSE
if("SMC_Contractile_Score" %in% colnames(stData@meta.data)) {
  stData$is_SMC <- stData$SMC_Contractile_Score > quantile(stData$SMC_Contractile_Score, 0.3, na.rm=TRUE)
}
if("SMC_Synthetic_Score" %in% colnames(stData@meta.data) && !stData$is_SMC) {
  stData$is_SMC <- stData$SMC_Synthetic_Score > quantile(stData$SMC_Synthetic_Score, 0.3, na.rm=TRUE)
}


smc_meta <- stData@meta.data[stData$is_SMC, ]
non_smc_meta <- stData@meta.data[!stData$is_SMC, ]

if(!is.null(top_var_gene) && nrow(smc_meta) > 10) {
  cor_smc_zn_ecm <- cor.test(smc_meta$Zinc_Score, smc_meta$Top_ECM_Gene, method = "spearman")
  cor_non_smc_zn_ecm <- cor.test(non_smc_meta$Zinc_Score, non_smc_meta$Top_ECM_Gene, method = "spearman")

  sink("ST10.Zn_ECM_Correlation.txt", append=TRUE)
  cat("
【SMC spot中锌-", top_var_gene, "相关性】
")
  cat("  SMC spot数:", nrow(smc_meta), "
")
  cat("  Spearman rho:", round(cor_smc_zn_ecm$estimate, 4), "
")
  cat("  P-value:", format(cor_smc_zn_ecm$p.value, digits=4, scientific=TRUE), "

")
  cat("【非SMC spot中锌-", top_var_gene, "相关性】
")
  cat("  非SMC spot数:", nrow(non_smc_meta), "
")
  cat("  Spearman rho:", round(cor_non_smc_zn_ecm$estimate, 4), "
")
  cat("  P-value:", format(cor_non_smc_zn_ecm$p.value, digits=4, scientific=TRUE), "
")
  sink()
}


if(!is.null(top_var_gene) && all(c("Zn-high", "Zn-low") %in% stData$Zn_Level)) {


  kw_ecm_module <- kruskal.test(ECM_Degradation_Score ~ Zn_Level, data = stData@meta.data)


  kw_ecm_top <- kruskal.test(Top_ECM_Gene ~ Zn_Level, data = stData@meta.data)


  p_ecm_module_box <- ggplot(stData@meta.data, aes(x = Zn_Level, y = ECM_Degradation_Score, fill = Zn_Level)) +
    geom_boxplot() +
    scale_fill_manual(values = c("Zn-high" = "#E74C3C", "Zn-mid" = "#95A5A6", "Zn-low" = "#3498DB")) +
    ggtitle(paste("ECM Module by Zn Level
K-W P =", format(kw_ecm_module$p.value, digits=3))) +
    theme_minimal() +
    theme(legend.position = "none")


  p_ecm_top_box <- ggplot(stData@meta.data, aes(x = Zn_Level, y = Top_ECM_Gene, fill = Zn_Level)) +
    geom_boxplot() +
    scale_fill_manual(values = c("Zn-high" = "#E74C3C", "Zn-mid" = "#95A5A6", "Zn-low" = "#3498DB")) +
    ggtitle(paste(top_var_gene, "by Zn Level
K-W P =", format(kw_ecm_top$p.value, digits=3))) +
    theme_minimal() +
    theme(legend.position = "none")

  pdf(file="ST10.Zn_ECM_BoxPlot.pdf", width=12, height=5)
  print(p_ecm_module_box + p_ecm_top_box)
  dev.off()


  sink("ST10.Zn_ECM_KWtest.txt")
  cat("=== 锌区域间ECM差异（Kruskal-Wallis） ===

")
  cat("【ECM模块评分】
")
  cat("  H =", round(kw_ecm_module$statistic, 4), ", P =", format(kw_ecm_module$p.value, digits=4), "

")
  cat("【", top_var_gene, "】
")
  cat("  H =", round(kw_ecm_top$statistic, 4), ", P =", format(kw_ecm_top$p.value, digits=4), "
")
  sink()
}


if(!is.null(top_var_gene)) {
  stData$Zn_ECM_Hotspot <- (stData$Zinc_Score >= quantile(stData$Zinc_Score, 0.75, na.rm=TRUE)) &
                            (stData$Top_ECM_Gene >= quantile(stData$Top_ECM_Gene, 0.75, na.rm=TRUE))

  p_hotspot <- SpatialDimPlot(stData, group.by = "Zn_ECM_Hotspot", pt.size.factor = 1.5) +
    scale_fill_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#E8E8E8"), name = "Zn-ECM Hotspot") +
    ggtitle(paste("Zn-", top_var_gene, "Hotspots"))

  pdf(file="ST10.ZnECM_Hotspot_Spatial.pdf", width=10, height=8)
  print(p_hotspot)
  dev.off()
}


if(!is.null(top_var_gene)) {
  p_scatter <- ggplot(stData@meta.data, aes(x = Zinc_Score, y = Top_ECM_Gene, color = Zn_Level)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    scale_color_manual(values = c("Zn-high" = "#E74C3C", "Zn-mid" = "#95A5A6", "Zn-low" = "#3498DB")) +
    ggtitle(paste("Zinc Score vs", top_var_gene, "
Spearman rho =", round(cor_zn_ecm_top$estimate, 3))) +
    theme_minimal()

  pdf(file="ST10.Zn_vs_ECM_Scatter.pdf", width=10, height=8)
  print(p_scatter)
  dev.off()
}

if(!is.null(top_var_gene)) {
  p_top_gene <- SpatialFeaturePlot(stData, features = "Top_ECM_Gene", pt.size.factor = 1.5) +
    scale_fill_gradientn(colors = c("darkblue", "blue", "lightblue", "white", "yellow", "orange", "red"),
                         name = top_var_gene) +
    ggtitle(paste(top_var_gene, "Spatial Distribution")) +
    theme(legend.position = "right")

  pdf(file=paste0("ST10.", top_var_gene, "_Spatial.pdf"), width=10, height=8)
  print(p_top_gene)
  dev.off()
}


save(stData, file="ST_with_Zn_ECM_Analysis_v2.Rdata")


sink("ST10.Zn_ECM_Summary_v2.txt")
cat("=== 锌离子空间差异与ECM降解关联分析总结（修正版）===

")
cat("样本:", sampleName, "
")
cat("总spot数:", ncol(stData), "
")
cat("锌相关基因检出:", length(zinc_genes_present), "个
")
cat("ECM相关基因检出:", length(ecm_present), "个
")
if(!is.null(top_var_gene)) {
  cat("Top变异ECM基因:", top_var_gene, "
")
}
cat("
")

cat("【锌离子区域分布】
")
print(table(stData$Zn_Level))
cat("
")

cat("【锌评分范围】
")
cat("  Min:", round(min(stData$Zinc_Score, na.rm=TRUE), 4), "
")
cat("  Max:", round(max(stData$Zinc_Score, na.rm=TRUE), 4), "
")
cat("  Mean:", round(mean(stData$Zinc_Score, na.rm=TRUE), 4), "

")

cat("【SMC spot定义】
")
cat("  SMC spot数:", sum(stData$is_SMC), "
")
cat("  占比:", round(sum(stData$is_SMC)/ncol(stData)*100, 2), "%

")

if(!is.null(top_var_gene)) {
  cat("【锌-Top基因相关性】
")
  cat("  全样本 Spearman rho:", round(cor_zn_ecm_top$estimate, 4), "
")
  cat("  P-value:", format(cor_zn_ecm_top$p.value, digits=4, scientific=TRUE), "

")

  cat("【锌-ECM热点】
")
  cat("  热点spot数:", sum(stData$Zn_ECM_Hotspot), "
")
  cat("  占比:", round(sum(stData$Zn_ECM_Hotspot)/ncol(stData)*100, 2), "%

")
}

sink()

