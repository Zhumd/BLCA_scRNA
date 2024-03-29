#step1 load data and packages
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
outdir = '/singlecell/CeleScope/LHBC_all/08.analysis/'
sample = 'LHBCall' 
load(paste0(outdir,"LHBC0911.RData"))
load(paste0(outdir,"LHBC0917.RData"))
load(paste0(outdir,"LHBC0803.RData"))
load(paste0(outdir,"LH0826.RData"))
rds.list <- list(LHBC0803.rds, LHBC0917.rds, LHBC0911.rds, LH0826.rds)
#normalize and identify variable features for each dataset independently
rds.list <- lapply(rds.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = 'vst', nfeatures=2000)
})
#select features that are repeatedly variable across datasets for integration.
features <- SelectIntegrationFeatures(object.list = rds.list)

#step2 identify anchors and use these anchors to integrate these datasets.
bladder.anchors <- FindIntegrationAnchors(rds.list, anchor.features = features)
bladder.combined <- IntegrateData(bladder.anchors)
DefaultAssay(bladder.combined) <- 'integrated'
#run the standard workflow for visualization and clustering.
bladder.combined <- ScaleData(bladder.combined, verbose = FALSE)
bladder.combined <- RunPCA(bladder.combined, npcs = 30, verbose = FALSE)
bladder.combined <- RunUMAP(bladder.combined, dims = 1:30, reduction = 'pca')
bladder.combined <- FindNeighbors(bladder.combined, reduction='pca', dims = 1:30)
bladder.combined <- FindClusters(bladder.combined, resolution = 0.5) 
DefaultAssay(bladder.combined) <- 'RNA'
pdf(paste0(outdir,'umap_', sample, '.pdf'),width = 15)
p1 <- DimPlot(bladder.combined, reduction = 'umap', group.by = 'orig.ident')
p2 <- DimPlot(bladder.combined, reduction = 'umap', label = TRUE,
              repel = TRUE)
p1+p2
dev.off()
pdf(paste0(outdir,'umap_split_', sample, '.pdf'), width = 15)
DimPlot(bladder.combined, reduction = 'umap', split.by = 'orig.ident')
dev.off()

#step3 DEG analysis, find markers for every cluster compared to all remaining cells, report only the positive ones
all.markers <- FindAllMarkers(bladder.combined, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25 )
write.csv(all.markers, paste0(outdir, 'allClusterMarkers_',sample, '.csv'))
top10 <- all.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
top10order1 <- top10 %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group=TRUE) #descending order by cluster.
top4marker <- top10order1 %>% group_by(cluster) %>% top_n(n=4, wt=avg_log2FC)
top4marker[1:9, c('cluster','gene')]
write.csv(top4marker, paste0(outdir, 'allClusterMarkers_top4_',sample, '.csv'))
           
markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5", 
                     "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1", 
                     "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")


markers.to.plot <- unique(markers.to.plot)
pdf(paste0(outdir,'dot_top4genes_', sample, '.pdf'), height  = 15)
DotPlot(bladder.combined, features=markers.to.plot, cols = c("white","blue"),
        dot.scale = 8)+ggplot2::coord_flip() 
dev.off()

genes_to_check = c("PTPRC","EPCAM","CD3G","CD3E", "CD79A", "BLNK","MS4A1", "CD68", "CSF1R", 
                   "MARCO", "CD207", "PMEL", "ALB", "C1QB", "CLDN5", "FCGR3B", "COL1A1")
pdf(paste0(outdir,'genes_to_check_', sample, '.pdf'), height  = 7)
genes_to_check = c("PTPRC","EPCAM","CD3G","CD3E", "CD79A", "BLNK","MS4A1", "CD68", "CSF1R", 
                   "MARCO", "CD207", "PMEL", "ALB", "C1QB", "CLDN5", "FCGR3B", "COL1A1")
pdf(paste0(outdir,'genes_to_check_', sample, '.pdf'), height  = 7)
# All on Dotplot 
DotPlot(bladder.combined, features = genes_to_check) +ggplot2::coord_flip()
dev.off()
