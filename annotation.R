library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
outdir = '/singlecell/CeleScope/LHBC_all/08.analysis/'
sample = 'LHBCall' 
DefaultAssay(bladder.combined) <- 'RNA'
epi <- c('0','1','2','3','4','5','6','7','11','14')
imm <- c('8','10','16','13')
stromal <- c('9','12','15')
bladder.combined@meta.data$immune_annotation <- ifelse(bladder.combined@meta.data$seurat_clusters %in% imm, 'immune',
                                                       ifelse(bladder.combined@meta.data$seurat_clusters %in% epi, 'epithelial', 'stromal'))

new <- paste(bladder.combined@meta.data$immune_annotation,
             bladder.combined@meta.data$seurat_clusters)
bladder.combined@meta.data$new <- new
genes_to_check=c('PTPRC','CD3G','CD3E','CD79A','BLNK','CD68','CSF1R','MARCO','CD207','PMEL','MLANA','PECAM1','CD34','VWF','EPCAM','SFN','KRT19','ACTA2','MCAM','MYLK','MYL9','FAP','THY1','ALB')
pdf(paste0(outdir,'general_cell_marker_', sample, '.pdf'),width = 15)
DotPlot(bladder.combined, features = genes_to_check,
        assay='RNA',group.by = 'new' )+RotatedAxis() 
dev.off()
Malepi <- infercnvEpitheial[infercnvEpitheial$infercnv_annotation=="tumor",]$cell_id
Macro <- c('8','10','16')
T <- c('13')
Fibroblasts <- c('9','15') 
Endothelial <- c('12')
bladder.combined@meta.data$immune_annotation2 <- ifelse(bladder.combined@meta.data$seurat_clusters %in% Macro, 'Macrophage',
                                                        ifelse(bladder.combined@meta.data$seurat_clusters %in% T, 'T cells',
                                                               ifelse(bladder.combined@meta.data$seurat_clusters %in% Fibroblasts, 'Fibroblasts',
                                                                      ifelse(bladder.combined@meta.data$seurat_clusters %in% Endothelial, 'Endothelial',
                                                                             ifelse(rownames(bladder.combined@meta.data) %in% Malepi, 'Maligant-epithelial','Nonmaligant-epithelial')))))
infercnvEpitheial <- infercnv
infercnvEpitheial$infercnv_annotation2 <- ifelse(infercnvEpitheial$infercnv_label %in% c(3,4),'nontumor', 'tumor')
Malepi <- infercnvEpitheial[infercnvEpitheial$infercnv_annotation2=="tumor",]$cell_id
Macro <- c('8','10','16')
T <- c('13')
Fibroblasts <- c('9','15') 
Endothelial <- c('12')
bladder.combined@meta.data$immune_anno_cluster34NT <- ifelse(bladder.combined@meta.data$seurat_clusters %in% Macro, 'Macrophage',
                                                             ifelse(bladder.combined@meta.data$seurat_clusters %in% T, 'T cells',
                                                                    ifelse(bladder.combined@meta.data$seurat_clusters %in% Fibroblasts, 'Fibroblasts',
                                                                           ifelse(bladder.combined@meta.data$seurat_clusters %in% Endothelial, 'Endothelial',
                                                                                  ifelse(rownames(bladder.combined@meta.data) %in% Malepi, 'Maligant-epithelial','Nonmaligant-epithelial')))))

# save data
save(bladder.combined,file= paste0(outdir,'thrid_bladderCombined.Rdata'))
phe=bladder.combined@meta.data
save(phe, file=paste0(outdir,'phe_of_immune_anno_cluster34NT.Rdata'))
bladder.combined@meta.data$immune_annotation3 <- ifelse(bladder.combined@meta.data$immune_annotation2=='Maligant-epithelial','ME',
                                                        ifelse(bladder.combined@meta.data$immune_annotation2=='Nonmaligant-epithelial','NME','Others'))
bladder.combined@meta.data$es <- colMeans(bladder.combined[rownames(bladder.combined) %in% genes_to_check, ])
pdf(paste0(outdir,'epithelial_score_', sample, '.pdf'))
VlnPlot(bladder.combined, features = 'es', pt.size = 0, group.by = 'immune_annotation3')+
  ylab('Epithelial score')+xlab('')+ggtitle('')+
  stat_summary(fun = median, geom = 'point', size = 10, colour = 'black', shape = 95)+NoLegend()
dev.off()
#heat map
general_markers <- data.frame('gene'=c('EPCAM','SFN','KRT19','PECAM1','CD34','VWF','ACTA2','MCAM','MYLK','MYL9','FAP','THY1','PTPRC','CD3G','CD3E','CD79A','BLNK','CD68','CSF1R','MARCO','CD207' ),
                              'type'=c(rep('Epithelial',3 ),rep('Endothelial',3 ),rep('Fibroblasts',6),rep('Immune',9)))

install.packages("broman")
library(broman)
epithelial_col <- brocolors("crayons")["Maroon"]
fibroblast_col <- brocolors("crayons")["Aquamarine"]
endothelial_col <- brocolors("crayons")["Wisteria"]
immune_col <- brocolors("crayons")["Inchworm"]
marker_cols <- c("epithelial" = unname(epithelial_col), 
                 "fibroblast" = unname(fibroblast_col), "endothelial" = unname(endothelial_col), 
                 "immune" = unname(immune_col)
)
cycling_mel_cols <- c("non-cycling" = "gainsboro",
                      "cycling" = unname(brocolors("crayons")["Mulberry"]))
pats_cols <- c("LH0826" = unname(brocolors("crayons")["Orange Red"]), "LHBC0803" = unname(brocolors("crayons")["Orange"]), 
                             "LHBC0911" = unname(brocolors("crayons")["Blue Violet"]), "LHBC0917" = unname(brocolors("crayons")["Sky Blue"]))
tsne_cols <- c("epithelial" = unname(basal_epithelial_col), "stroma" = unname(stroma_col), "endothelial" = unname(endothelial_col),
               "Tcell" = unname(t_cell_col), "Bcell" = unname(b_cell_col), "macrophage" = unname(macrophage_col))
anno_colors <- list("marker" = marker_cols, "cycling" = cycling_mel_cols, 
                    "patient" = pats_cols)

general_markers <-data.frame('gene'=c('EPCAM','SFN','KRT19','PECAM1','CD34','VWF','ACTA2','MCAM','MYLK','MYL9','FAP','THY1','PTPRC','CD3G','CD3E','CD79A','BLNK','CD68','CSF1R','MARCO','CD207' ),
           'type'=c(rep('epithelial',3 ),rep('endothelial',3 ),rep('fibroblast',6),rep('immune',9)), stringsAsFactors = FALSE)

markers <- general_markers
markers$type_heatmap <- markers$type
colors_markers_ch <- markers$type_heatmap
for (i in c(1:length(names(anno_colors$marker)))) {
  colors_markers_ch <- replace(colors_markers_ch, colors_markers_ch == names(anno_colors$marker)[i], anno_colors$marker[i])
  }
splits_ch <- as.factor(markers$type_heatmap) #明白了，这是为了row_title，正好是四个大类
splits_ch <- factor(splits_ch, levels(splits_ch)[c(2,3,1,4)])
head(splits_ch)
colors_anno_markers_ch <- splits_ch
ha_rows <- HeatmapAnnotation(df = data.frame(type=colors_anno_markers_ch),
                             annotation_legend_param = list(type = list(ncol = 2, title = "cell type", title_position = "topcenter")),
                             which = "row", col = list("type"=marker_cols),annotation_width = unit(3, "mm"))
mats_nowall <- as.matrix(bladder.combined@assays$integrated@data)
mats_now <- mats_nowall[,sample(1:ncol(mats_now), 2000, replace = FALSE)]
dim(mats_now)
a=ha_rows + 
  Heatmap(mats_now[match(as.character(markers$gene), rownames(mats_now)),], 
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE,
          name = 'patients_now[1]', clustering_distance_columns = "euclidean", row_names_side = "right", 
          row_names_gp = gpar(fontsize = 10, col = colors_markers_ch),
          split = splits_ch, gap = unit(1, "mm"), column_title = 'patients_now[1]', 
          column_title_gp = gpar(fontsize = 11),
          row_title_gp = gpar(font = 11),
          heatmap_legend_param = list(title = "expression", title_position = "topcenter", color_bar = "continuous"),
 )
pdf(paste0(outdir,'markers_heatmap_', sample, '.pdf'), height  = 7,width = 8)
draw(a, gap = unit(0.1, "cm"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
#annotation bar
mats_now <- mats_nowall[,sample(1:ncol(mats_now), 2000, replace = FALSE)]
phe_now <- phe[rownames(phe) %in% colnames(mats_now),]
patient <-data.frame(patient=phe_now$orig.ident,row.names=rownames(phe_now))
#make cycling and non-cycling data 
melanoma_cellcycle <- read.table(paste0(outdir, "melanoma_cellcycle.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(melanoma_cellcycle) <- c('G1S','G2M')
melanoma_g1s <- melanoma_cellcycle$G1S
# function to average the expression over a vector of gene indices, for all cells
avg_expr_genes <- function(mat_to_fit, indices){
  scores_cells <- apply(mat_to_fit, 2, function(x){mean(x[indices])})
  return(scores_cells)
}
# function to match and clean a vector of genes
match_clean_vector_genes <- function(mat_to_fit, vec_genes){
  
  vec_genes <- unique(vec_genes)
  vec_genes <- trimws(vec_genes)
  vec_genes <- data.frame("gene" = vec_genes)
  vec_genes$index <- match(vec_genes$gene, rownames(mat_to_fit))
  if (length(which(is.na(vec_genes$index)) > 0))
    vec_genes <- vec_genes[-which(is.na(vec_genes$index)), ]
  return(vec_genes)
}
mats_nowall <- as.matrix(bladder.combined@assays$integrated@data)
melanoma_g1s <- match_clean_vector_genes(mats_nowall, melanoma_g1s)
scores_g1s <- avg_expr_genes(mats_nowall, melanoma_g1s$index)
melanoma_g2m <- melanoma_cellcycle$G2M
melanoma_g2m <- match_clean_vector_genes(mats_nowall, melanoma_g2m)
scores_g2m <- avg_expr_genes(mats_nowall, melanoma_g2m$index)
cycling_mel <- rep(NA, length(scores_g1s))
for (i in 1:length(cycling_mel)) {
  if (scores_g1s[i] >= (median(scores_g1s) + 2 * mad(scores_g1s)) && scores_g2m[i] < (median(scores_g2m) + 2 * mad(scores_g2m)))
    cycling_mel[i] <- "cycling"
  if (scores_g1s[i] < (median(scores_g1s) + 2 * mad(scores_g1s)) && scores_g2m[i] >= (median(scores_g2m) + 2 * mad(scores_g2m)))
    cycling_mel[i] <- "cycling"
  if (scores_g1s[i] < (median(scores_g1s) + 2 * mad(scores_g1s)) && scores_g2m[i] < (median(scores_g2m) + 2 * mad(scores_g2m)))
    cycling_mel[i] <- "non-cycling"
  if (scores_g1s[i] >= (median(scores_g1s) + 2 * mad(scores_g1s)) && scores_g2m[i] >= (median(scores_g2m) + 2 * mad(scores_g2m)))
    cycling_mel[i] <- "cycling"
}
phe$mel_scores_g1s <- scores_g1s
phe$mel_scores_g2m <- scores_g2m
phe$cycling_mel <- cycling_mel
phe_now <- phe[rownames(phe) %in% colnames(mats_now),]
#save data
save(phe, file=paste0(outdir,'phe_of_immune_anno_cluster34NT_cycling.Rdata')) 
ha_cols_up <- HeatmapAnnotation(df=data.frame(cycling=phe_now$cycling_mel),
                                   col=list(cycling=anno_colors$cycling),
                                   which='column',
                                   show_annotation_name=TRUE,
                                   annotation_name_side='right',
                                   annotation_name_gp=gpar(fontsize=11),
                                   annotation_legend_param=list(list(title_position="topcenter",
                                                                     title=c('patient'))),
                                annotation_name_align = FALSE,
                                gap = unit(c(1,1), 'mm')
                                )
markers <-unique(general_markers[general_markers$gene %in% rownames(mats_now),])
markers$type_heatmap <- markers$type
colors_markers_ch <- markers$type_heatmap
colors_markers_ch <- as.character(colors_markers_ch)
for (i in c(1:length(names(anno_colors$marker)))) {
  colors_markers_ch <- replace(colors_markers_ch, colors_markers_ch == names(anno_colors$marker)[i], anno_colors$marker[i])
}
splits_ch <- as.factor(markers$type_heatmap) 
splits_ch <- factor(splits_ch, levels(splits_ch)[c(2,3,1,4)])
head(splits_ch)
colors_anno_markers_ch <- splits_ch
ha_rows <- HeatmapAnnotation(df = data.frame(type=colors_anno_markers_ch),
                             annotation_legend_param = list(type = list(ncol = 1, title = "cell type", title_position = "topcenter")),
                             which = "row", col = list("type"=marker_cols),annotation_width = unit(3, "mm"))

a=ha_rows + 
  Heatmap(mats_now[match(as.character(markers$gene), rownames(mats_now)),], 
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE,
          name = 'patients_now[1]', clustering_distance_columns = "euclidean", row_names_side = "right", 
          row_names_gp = gpar(fontsize = 10, col = colors_markers_ch),
          split = splits_ch, gap = unit(1, "mm"), column_title = 'patients_now[1]', 
          column_title_gp = gpar(fontsize = 11),
          row_title_gp = gpar(font = 11), top_annotation = ha_cols_up,
          heatmap_legend_param = list(title = "expression", title_position = "topcenter", color_bar = "continuous"),
 )
pdf(paste0(outdir,'markers_heatmap4_', sample, '.pdf'), height  = 7,width = 8)
draw(a, gap = unit(0.1, "cm"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

