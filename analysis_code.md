# Smart-CPSP Analysis Pipeline

## 1. SingleR Annotation (reference-based)

```r
library(celldex)
library(SingleR)
library(tidyverse)
ref <- fetchReference("mouse_rnaseq", "2024-02-26")

scRNA=readRDS("scRNA.v1.rds")
DefaultAssay(scRNA)<-'RNA'
matirx=GetAssayData(scRNA)
pred.scRNA <- SingleR(test = matirx, ref = ref, assay.type.test=1,
    labels = ref$label.main)
singleR_output=as.data.frame(pred.scRNA)
colnames(singleR_output)
singleR_output$label.score=str_replace_all(paste("scores",singleR_output$pruned.labels),' ','.')
df=singleR_output
df$scores.label <- df[cbind(1:nrow(df), match(df$label.score, names(df)))]
out=df[colnames(scRNA),c("pruned.labels",'scores.label')]
colnames(out)=c("singleR.labels",'singleR.score')
write.csv(out,"singleR1.cell.results.csv")

ref <- fetchReference("immgen", "2024-02-26")
pred.scRNA <- SingleR(test = matirx, ref = ref, assay.type.test=1,
    labels = ref$label.main)
singleR_output=as.data.frame(pred.scRNA)
singleR_output$label.score=str_replace_all(paste("scores",singleR_output$pruned.labels),' ','.')
df=singleR_output
colnames(df)
df$scores.label <- df[cbind(1:nrow(df), match(df$label.score, names(df)))]
out=df[colnames(scRNA),c("pruned.labels",'scores.label')]
colnames(out)=c("singleR.labels",'singleR.score')
write.csv(out,"singleR2.cell.results.csv")
```

---

## 2. Cell Clustering (SCTransform + PCA + UMAP + resolution)

```r
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
library(tidyverse)

obj=readRDS("scRNA.v1.rds")
obj$Sample=obj$Group
obj@meta.data<-obj@meta.data[,c("orig.ident",'nCount_RNA','nFeature_RNA','Group','Sample','percent.mt','percent.rpls','percent.hbb')]
obj <- SCTransform(obj, vst.flavor = "v2")
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 0.2)
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca",verbose =F)
singleR_anno=read.csv("singleR1.cell.results.csv",row.names = 1)
singleR_anno$singleR.labels=ifelse(singleR_anno$singleR.score>0.3,singleR_anno$singleR.labels,NA)
labels.counts=table(singleR_anno$singleR.labels)
rare_labels <- names(labels.counts[labels.counts < 10])
singleR_anno$singleR.labels[singleR_anno$singleR.labels %in% rare_labels] <- NA
obj=AddMetaData(obj,singleR_anno)
singleR_anno2=read.csv("singleR2.cell.results.csv",row.names = 1)
singleR_anno2$singleR.labels=ifelse(singleR_anno2$singleR.score>0.3,singleR_anno2$singleR.labels,NA)
labels.counts=table(singleR_anno2$singleR.labels)
rare_labels <- names(labels.counts[labels.counts < 10])
singleR_anno2$singleR.labels[singleR_anno2$singleR.labels %in% rare_labels] <- NA
colnames(singleR_anno2)=c("immgen.singleR.labels",'immgen.singleR.score')
obj=AddMetaData(obj,singleR_anno2)

#obj$predicted.subclass=ifelse(obj$predicted.subclass.score>0.3,obj$predicted.subclass,NA)
obj <- FindClusters(obj, resolution = c(0.4,0.5))
p1<-DimPlot(obj, reduction = 'umap', group.by = "SCT_snn_res.0.4",label=T,label.size = 6)
p2<-DimPlot(obj, reduction = 'umap', group.by = "SCT_snn_res.0.5",label=T,label.size = 6)
wrap_plots(p1, p2,widths = c(1,1))
p1<-DimPlot(obj, reduction = 'umap', group.by = "Group",label=F,label.size = 6)
p2<-DimPlot(obj, reduction = 'umap', group.by = "seurat_clusters",label=T,label.size = 6)
wrap_plots(p1, p2,widths = c(1,1))
p1<-DimPlot(obj, reduction = 'umap', group.by = "seurat_clusters",label=T,label.size = 6)
p2<-DimPlot(obj, reduction = 'umap', group.by = "singleR.labels",label=T,label.size = 5)
p3<-DimPlot(obj, reduction = 'umap', group.by = "immgen.singleR.labels",label=T,label.size = 5)
p1<-DimPlot(obj, reduction = 'umap', group.by = "seurat_clusters",label=T,label.size = 4,split.by="Sample",ncol=4)

library(scales)
QCplot(obj,group.by="seurat_clusters")
plot_list <- lapply(levels(obj), function(cluster) {
    cells_to_highlight <- Cells(obj)[obj$seurat_clusters == cluster]
    identities <- levels(Idents(obj))
    my_color_palette <- hue_pal()(length(identities))
    names(my_color_palette)=levels(Idents(obj))
    col=my_color_palette[cluster]
    DimPlot(obj, cells.highlight = cells_to_highlight,cols.highlight=col, reduction = 'umap',group.by="seurat_clusters") +
    ggtitle(paste("Cluster", cluster)) +
    theme(legend.position = "none")
})
combined_plot <- wrap_plots(plot_list, ncol = 7)
p1<-DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",label=T,label.size = 6)
p2<-DimPlot(obj, reduction = "umap", group.by = "singleR.labels",label=T,label.size = 6)
wrap_plots(p1, p2,widths = c(1,1))

# cell type
markers<-c(
    "Pecam1",'Cldn5','Cdh5',"Itm2a","Ly6c1","Pltp",#Endothelial
    "Kcnj8","Cald1", "Vtn", "Notch3",
    'Col1a1','Dcn','Ttr',#Fibroblast/VLMC
    'Gfap','Aqp4','Aldh1l1',"Agt","Slc1a2","Gja1","Plpp3",#Astrocyte
    'Mbp','Mobp','Bcas1','C1ql1','Enpp2',"Plp1","Mog","Trf",#Oligodendrocytes
    'Olig1','Pdgfra',#OPCs
    'Cx3cr1','P2ry12','Tmem119',"Hexb","Siglech",#Microglia
    'Mrc1',#Macrophage
    "Tmem212","Ccdc153","Dynlrb2", "Rsph4a",
    "Pdgfra",#NG2
    "Col23a1",#tanycytes
    'Map2','Rbfox3',"Celf4","Snap25","Meg3", "Snhg11","Ndrg4",#Neurons
    'Slc17a6','Slc17a7',#Excitatory neurons
    "Slc32a1","Gad2",'Gad1',#Inhibitory neurons
    "Mrc1","Cx3cr1",#ImmuneCells
    "Flt1",#Vascula
    "Notch1","Sox2",#stem cell
    "Epcam"#Epithelial cell
)

ave_exp = AverageExpression(object = obj2, assays = "RNA", features = markers)
markers=rownames(ave_exp$RNA[rowSums(ave_exp$RNA)>1,])
FeaturePlot(obj, features = markers ,label=T,ncol=7,reduction = 'umap')

clustermarkers <- FindAllMarkers(object = obj, only.pos = T,assay = 'SCT')
write.csv(clustermarkers,'main.cluster.markers.csv')

clustermarkers=read.csv("main.cluster.markers.csv",row.names = 1)
cDEG=clustermarkers%>%
    filter(p_val_adj<0.001,avg_log2FC>1,pct.1>0.5,pct.2<0.5)%>%
    arrange(cluster,desc(avg_log2FC))%>%
    group_by(cluster)%>%
    top_n(n=20,wt=avg_log2FC)
p<-DoHeatmap(object = subset(obj, downsample = 1000),features = cDEG$gene)

# Split the dataframe by cluster
split_df <- split(cDEG, cDEG$cluster)

# Extract the 'gene' column from each subset
gene_lists <- lapply(split_df, function(subset) {
  return(subset$gene)
})
print(gene_lists)

saveRDS(obj,"scRNA.v2.rds")
```

---

## 3. Main Cell Type Annotation

```r
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
library(tidyverse)

obj=readRDS("scRNA.v2.rds")
celltype=data.frame(ClusterID=0:19,
                  celltype= 0:19)
celltype[celltype$ClusterID %in% c(3,0),2]='Microglia'
celltype[celltype$ClusterID %in% c(6,10),2]='Astrocyte'
celltype[celltype$ClusterID %in% c(2,16),2]='OPCs'
celltype[celltype$ClusterID %in% c(8),2]='OL-diff'
celltype[celltype$ClusterID %in% c(1,4,5),2]='OL-mature'
celltype[celltype$ClusterID %in% c(19),2]='ECs'
celltype[celltype$ClusterID %in% c(17),2]='T cells'
celltype[celltype$ClusterID %in% c(9),2]='Macrophages'
celltype[celltype$ClusterID %in% c(18,15),2]='DCs'
celltype[celltype$ClusterID %in% c(7,11,13,14),2]='Neuron'
celltype[celltype$ClusterID %in% c(12),2]='Pericyte'
obj=celltype_rename(obj,celltype)

saveRDS(obj,"scRNA.anno.rds")
p1 <- DimPlot(obj, reduction = "umap", group.by='celltype',label=T,label.size = 5)#+xlim(-15,20)
p2 <- DimPlot(obj, reduction = "umap", group.by='seurat_clusters',label=T)#+xlim(-15,20)
old_info=readRDS("../analysis_batch05/scRNA.v2.rds")
old_batch=unique(old_info$Sample)
obj=subset(obj,Sample%!in%old_batch & celltype!="Neuron",invert=T)
obj$Group=obj$Sample
obj$Group=ifelse(obj$Group%in%c('saline','WT-saline'),'saline',obj$Group)
obj$Group=ifelse(obj$Group%in%c("SPP1-ko-control",'Spp1-ko-saline'),'SPP1-ko-saline',obj$Group)
obj$Group=ifelse(obj$Group%in%c("SPP1-ko-cpsp",'Spp1-ko-TH'),'SPP1-ko-TH',obj$Group)
obj$Group=ifelse(obj$Group%in%c("CPSP-day7",'WT-TH'),'WT-TH',obj$Group)
p1 <- DimPlot(obj, reduction = "umap", group.by='celltype',label=T,label.size = 5)+NoLegend()
ggsave("figs/main.celltype.pdf")

p1 <- DimPlot(obj, reduction = "umap", group.by='celltype',label=F,label.size = 3,split.by="Group",ncol=4)#+xlim(-15,20)
ggsave("figs/Sample.celltype.umap.pdf",width=16,height=8)

library(RColorBrewer)
library(scales)

plot_gene_optimized <- function(s.in, gene, min_cutoff=0,max_cutoff=NULL,
                                     ncol = 2, bg_color = "#F8F8F8",
                                     zero_color = "#D3D3D3") {
  fpm_matrix=GetAssayData(s.in,layer = "data",assay='SCT')
  tsne_coords=s.in@reductions$umap@cell.embeddings
  base_df <- data.frame(
    tSNE1 = tsne_coords[, 1],
    tSNE2 = tsne_coords[, 2]
  )
    expr <- as.numeric(fpm_matrix[gene, ])
    zero_idx <- expr < min_cutoff
    pos_idx <- expr >= min_cutoff

    if (!is.null(max_cutoff)) {
        expr[expr>max_cutoff]=max_cutoff
      }
    plot_data <- data.frame(
      tSNE1 = c(base_df$tSNE1[zero_idx], base_df$tSNE1[pos_idx]),
      tSNE2 = c(base_df$tSNE2[zero_idx], base_df$tSNE2[pos_idx]),
      Group = c(rep("Zero", sum(zero_idx)), rep("Positive", sum(pos_idx))),
      Expr = c(rep(0, sum(zero_idx)), expr[pos_idx] )
    )
    expr_rate <- round(mean(expr >= min_cutoff) * 100, 1)
    p <- ggplot() +
      geom_point(data = subset(plot_data, Group == "Zero"),
                 aes(x = tSNE1, y = tSNE2),
                 color = zero_color,
                 size = 0.3,
                 alpha = 0.4) +
      geom_point(data = subset(plot_data, Group == "Positive"),
                 aes(x = tSNE1, y = tSNE2, color = Expr),
                 size = 0.3,
                 alpha = 0.8) +
      scale_color_gradientn(
        colors = c("#FFFACD", "#FFD700", "#FFA500", "#FF6347", "#DC143C"),
        values = rescale(c(0, 1, 2, 3, 4)),
        name = "Expr",
        limits = c(min_cutoff, max(plot_data$Expr[plot_data$Group == "Positive"])),
        breaks = pretty_breaks(n = 4)
      ) +
      labs(
        title = gene,#paste0(gene, " (", expr_rate, "%)"),
        x = element_blank(),
        y = element_blank()
      ) +
      theme_minimal(base_size = 10) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        #plot.background = element_rect(fill = bg_color, color = NA),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        legend.key.height = unit(0.6, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        plot.margin = margin(2, 2, 2, 2, "mm"),
        aspect.ratio = 1
      )
}
p1<- plot_gene_optimized(obj,gene="WPRE",min_cutoff=4,max_cutoff=5)
p2<- plot_gene_optimized(obj,gene="Spp1",min_cutoff=4,max_cutoff=8)
ggsave(filename = "WPRE.featureplot.pdf",plot = p1,width = 6,height = 6)
ggsave(filename = "figs/Spp1.featureplot.pdf",plot = p2,width = 6,height = 6)

p1<- plot_gene_optimized(obj,gene="WPRE",min_cutoff=0.1,max_cutoff=5)
p2<- plot_gene_optimized(obj,gene="Spp1",min_cutoff=0.1,max_cutoff=8)

plot_proportion_barplot <- function(prop_df,celltype_threshold,title) {
    show_celltype=unique(prop_df$CellType[prop_df$sam_type_Proportion>celltype_threshold])
    plot_data=prop_df|>
        mutate(CellType=ifelse(CellType%in%show_celltype,prop_df$CellType,'other'))|>
        arrange(desc(sam_Proportion))
    identities <- unique(plot_data$CellType)
    colors <- rev(viridis::viridis(length(identities)))

    #colors <- hue_pal()(length(identities))
    plot_data$CellType=factor(plot_data$CellType,levels=identities)
    p <- ggplot(plot_data, aes(x = Sample, y = sam_Proportion, fill = CellType)) +
        geom_bar(stat = "identity", position = "stack", width = 0.7) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
              axis.text.y = element_text(size = 10),
              axis.title = element_text(size = 12, face = "bold"),
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              legend.position = "right",
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 9),
              aspect.ratio = 0.8
             ) +
          labs(
            title =title,
            x = element_blank(),
            y = "Proportion(%)"
          ) +
        scale_fill_manual(values = colors) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
    return(p)
}
get_positive_cell_prop <- function(seurat_obj,
                                   gene_name,
                                   sample_col="Sample",
                                   celltype_col = "celltype",
                                   assay = "RNA",
                                   slot = "data",
                                   threshold = 0,
                                   return_plot = FALSE) {
    expr_matrix <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
    gene_expr <- expr_matrix[gene_name, ]
    data <- data.frame(
        Sample = seurat_obj[[sample_col]][,1],
        CellType = seurat_obj[[celltype_col]][,1],
        Expression = gene_expr,
        stringsAsFactors = FALSE
    )
    is_positive <- gene_expr > threshold
    celltypes <- seurat_obj@meta.data[[celltype_col]]
    result <- data.frame()
    unique_celltypes <- unique(celltypes)
      sample_total_cells <- data %>%
        group_by(Sample) %>%
        summarise(Sample_Total_Cells = n(), .groups = 'drop')
    result <- data %>%
        mutate(is_positive = Expression > threshold) %>%
        group_by(Sample, CellType) %>%
        summarise(
            Positive_Count = sum(is_positive, na.rm = TRUE),
              sample_celltype_Cells = n(),
              sam_type_Proportion = Positive_Count / sample_celltype_Cells*100,
              #Positive_Percentage = Positive_Proportion * 100,
              .groups = 'drop'
        )|>
        left_join(sample_total_cells, by = "Sample") %>%
        mutate(sam_Proportion = Positive_Count / Sample_Total_Cells*100)
  return(result)
}
res <- get_positive_cell_prop(obj,sample_col="Group",celltype_col = "celltype",
                              assay = "SCT",slot = "data",
                              gene_name = "Spp1", threshold = 4,return_plot=F)
res2=res|>
    filter(!str_detect(Sample,'ko|cOpn'))
res2$Sample=factor(res2$Sample,levels=c("saline",'WT-TH','CPSP-3w'))
p<-plot_proportion_barplot(res2,celltype_threshold=10,title="Spp1")
ggsave(filename = "figs/Spp1.Proportion.pdf",plot = p,width = 6,height = 6)

sub_obj=subset(obj,Group%in%c("cOpn5-day1",'Cx3crl-cOpn5-3w','Cx3crl-cOpn5-No-Tamox'))
p1<- plot_gene_optimized(sub_obj,gene="WPRE",min_cutoff=1,max_cutoff=5)+ggtitle("min_cutoff 1")
p2<- plot_gene_optimized(sub_obj,gene="WPRE",min_cutoff=2,max_cutoff=5)+ggtitle("min_cutoff 2")
p1<- plot_gene_optimized(sub_obj,gene="WPRE",min_cutoff=3,max_cutoff=5)+ggtitle("min_cutoff 3")
p2<- plot_gene_optimized(sub_obj,gene="WPRE",min_cutoff=4,max_cutoff=5)+ggtitle("min_cutoff 4")
p1<- plot_gene_optimized(sub_obj,gene="WPRE",min_cutoff=3,max_cutoff=5)+ggtitle("WPRE(cOpn5)")
ggsave(filename = "figs/WPRE.featureplot.pdf",plot = p1,width = 6,height = 6)

SCT_matrix=GetAssayData(sub_obj,layer = "data",assay='SCT')
sub_obj@meta.data$WPRE_data=SCT_matrix["WPRE",rownames(sub_obj@meta.data)]
sub_obj@meta.data$WPRE_POS=ifelse(sub_obj@meta.data$WPRE_data>3,"pos",'neg')
data=as.data.frame.matrix(table(sub_obj@meta.data[,c("WPRE_POS",'celltype')]))

library(reshape2)
library(dplyr)
data <- data.frame(
  CellType = c("Astrocyte", "DCs", "ECs", "Macrophages", "Microglia",
               "Neuron", "OL-diff", "OL-mature", "OPCs", "Pericyte", "T cells"),
  neg = c(904, 15, 80, 77, 332, 1020, 69, 4027, 620, 200, 20),
  pos = c(55, 0, 15, 17, 1174, 4, 2, 34, 41, 11, 2)
)
data_long <- melt(data, id.vars = "CellType",
                  variable.name = "Status",
                  value.name = "Count")
data_long <- data_long %>%
  group_by(CellType) %>%
  mutate(Percentage = Count / sum(Count) * 100)
data_long$CellType=factor(data_long$CellType,levels=c("Microglia",'Macrophages','Astrocyte','OPCs','OL-mature','OL-diff','Neuron','Pericyte','ECs','T cells','DCs'))
p<-ggplot(data_long, aes(x = CellType, y = Percentage, fill = Status)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("neg" = "grey", "pos" = "#E64B35FF"),
                    labels = c("neg" = "Negative", "pos" = "Positive")) +
  labs(x = "Cell Type",
       y = "Percentage (%)",
       title = "Distribution of WPRE+ Cells Across Cell Types",
       fill = "Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "top")
ggsave(filename = "figs/WPRE+.celltype.percent.pdf",plot = p,width = 6,height = 5)
p <- ggplot(data_long, aes(x = CellType, y = Percentage, fill = Status)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(
    data = subset(data_long, Status == "pos"),
    aes(label = sprintf("%.1f", Percentage),
        y = Percentage + 2),
    size = 3.5,
    color = "black"
  ) +
  scale_fill_manual(values = c("neg" = "grey", "pos" = "#E64B35FF"),
                    labels = c("neg" = "Negative", "pos" = "Positive")) +
  labs(x = "Cell Type",
       y = "Percentage (%)",
       title = "Distribution of WPRE+ Cells Across Cell Types",
       fill = "Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "top")
ggsave(filename = "figs/WPRE+.celltype.percent.pdf",plot = p,width = 6,height = 5)
p<-ggplot(data_long, aes(x = CellType, y = Count, fill = Status)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("neg" = "grey", "pos" = "#E64B35FF"),
                    labels = c("neg" = "Negative", "pos" = "Positive")) +
  labs(x = "Cell Type",
       y = "Count",
       title = "Distribution of WPRE+ Cells Across Cell Types",
       fill = "Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "top")
ggsave(filename = "figs/WPRE+.celltype.num.pdf",plot = p,width = 6,height = 5)
percent_data <- data_long %>%
  group_by(CellType) %>%
  summarise(
    Total = sum(Count),
    Pos_Count = sum(Count[Status == "pos"]),
    Label = paste0(Pos_Count, "/", Total)
  )
p <- ggplot(data_long, aes(x = CellType, y = Count, fill = Status)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(
    data = percent_data,
    aes(x = CellType,
        y = Total + max(Total) * 0.05,
        label = Label),
    inherit.aes = FALSE,
    size = 3,
    color = "black"
  ) +
  scale_fill_manual(values = c("neg" = "grey", "pos" = "#E64B35FF"),
                    labels = c("neg" = "Negative", "pos" = "Positive")) +
  labs(x = "Cell Type",
       y = "Count",
       title = "Distribution of WPRE+ Cells Across Cell Types",
       fill = "Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "top")
ggsave(filename = "figs/WPRE+.celltype.num.pdf",plot = p,width = 6,height = 5)
pos_cells <- sub_obj$WPRE_POS == "pos"
celltype_prop <- prop.table(table(sub_obj$celltype[pos_cells]))
print(celltype_prop)
celltype_counts <- table(sub_obj$celltype[pos_cells])
celltype_prop <- prop.table(celltype_counts)
result <- data.frame(
  celltype = names(celltype_counts),
  count = as.vector(celltype_counts),
  proportion = as.vector(celltype_prop),
  percentage=  as.vector(celltype_prop)*100
)
print(result)
n_colors <- nrow(result)
if (n_colors <= 8) {
  colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
              "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
  colors <- colors[1:n_colors]
} else {
  colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_colors)
  # colors <- rainbow(n_colors)
}
ggplot(result, aes(x = 2, y = proportion, fill = celltype)) +
  geom_col(width = 1, color = "white", size = 0.8, alpha = 0.9) +
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(label = percentage),
            position = position_stack(vjust = 0.5),
            size = 4.5,
            fontface = "bold",
            color = "white",
            stroke = 0.2) +
  scale_fill_manual(values = colors) +
  xlim(0.5, 2.5) +
  theme_void() +
  labs(title = "Cell Type Distribution in WPRE_POS Positive Cells",
       subtitle = paste0("Total cells: ", sum(result$n), " (WPRE_POS positive)")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40", margin = margin(b = 15)),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  ) +
  guides(fill = guide_legend(reverse = FALSE, ncol = 1, byrow = TRUE))
result|>
    arrange(desc(proportion))
threshold=0.02
result_merge=result %>%
  mutate(celltype = ifelse(proportion < threshold & celltype!="Neuron", "other", celltype)) %>%
  group_by(celltype) %>%
  summarise(
    count = sum(count),
    proportion = sum(proportion),
    percentage = sum(percentage),
    .groups = "drop"
  ) %>%
  arrange(celltype != "other", celltype)
sum(as.numeric(data$percentage))
data <-result_merge
data <- data %>%
  arrange(desc(proportion)) %>%
  mutate(
    end_angle = 2 * pi * cumsum(proportion),
    start_angle = lag(end_angle, default = 0),
    mid_angle = (start_angle + end_angle) / 2,
    percent_label = sprintf("%.2f%%", percentage),
    celltype = factor(celltype, levels = celltype)
  )
colors <- c(
  "Microglia"   = "#E64B35FF",
  "Astrocyte"   = "#4DBBD5FF",
  "OPCs"        = "#00A087FF",
  "other"       = "#3C5488FF",
  "OL-mature"   = "#F39B7FFF",
  "Neuron"   = "#911EB4FF"
)
ggplot(data, aes(x = 2, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", size = 0.6) +
  coord_polar(theta = "y", start = 0) +

  #geom_text(aes(x = 2.3,
  #              label = percent_label,
  #              angle = 90 - 180 * mid_angle / pi),
  #          size = 4.5,
  #          fontface = "bold",
  #          color = "black") +
  scale_fill_manual(
    name = "",
    values = colors[levels(data$celltype)],
    breaks = levels(data$celltype),
    labels = paste0(levels(data$celltype), " (",
                    sprintf("%.2f%%", data$percentage[match(levels(data$celltype), data$celltype)]), ")")
  ) +
  xlim(0.5, 2.5) +
  labs(
    title = "WPRE_POS Positive Cells",
    subtitle = paste0("Total cells: ", sum(data$count))
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40", margin = margin(b = 15)),
    legend.position = "right",
    legend.text = element_text(size = 10.5),
    legend.key.size = unit(0.7, "cm"),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(filename = "figs/WPRE.donut_chart.0730.pdf",width = 8,height = 6)
SCT_matrix=GetAssayData(obj,layer = "data",assay='SCT')
obj@meta.data$Spp1_data=SCT_matrix["Spp1",rownames(obj@meta.data)]
obj@meta.data$Spp1_POS=ifelse(obj@meta.data$Spp1_data>3,"pos",'neg')
pos_cells <- obj$Spp1_POS == "pos"
celltype_prop <- prop.table(table(obj$celltype[pos_cells]))
print(celltype_prop)
celltype_counts <- table(obj$celltype[pos_cells])
celltype_prop <- prop.table(celltype_counts)
result <- data.frame(
  celltype = names(celltype_counts),
  count = as.vector(celltype_counts),
  proportion = as.vector(celltype_prop),
  percentage=  as.vector(celltype_prop)*100
)
print(result)
n_colors <- nrow(result)
if (n_colors <= 8) {
  colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
              "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
  colors <- colors[1:n_colors]
} else {
  colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_colors)
  # colors <- rainbow(n_colors)
}
ggplot(result, aes(x = 2, y = proportion, fill = celltype)) +
  geom_col(width = 1, color = "white", size = 0.8, alpha = 0.9) +
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(label = percentage),
            position = position_stack(vjust = 0.5),
            size = 4.5,
            fontface = "bold",
            color = "white",
            stroke = 0.2) +
  scale_fill_manual(values = colors) +
  xlim(0.5, 2.5) +
  theme_void() +
  labs(title = "SPP1_POS Positive Cells",
       subtitle = paste0("Total cells: ", sum(result$n), " (WPRE_POS positive)")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40", margin = margin(b = 15)),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  ) +
  guides(fill = guide_legend(reverse = FALSE, ncol = 1, byrow = TRUE))
threshold=0.02
result_merge=result %>%
  mutate(celltype = ifelse(proportion < threshold, "other", celltype)) %>%
  group_by(celltype) %>%
  summarise(
    count = sum(count),
    proportion = sum(proportion),
    percentage = sum(percentage),
    .groups = "drop"
  ) %>%
  arrange(celltype != "other", celltype)
data <- result_merge %>%
  arrange(desc(proportion)) %>%
  mutate(
    end_angle = 2 * pi * cumsum(proportion),
    start_angle = lag(end_angle, default = 0),
    mid_angle = (start_angle + end_angle) / 2,
    percent_label = sprintf("%.2f%%", percentage),
    celltype = factor(celltype, levels = celltype)
  )
colors <- c(
  "Microglia"   = "#E64B35FF",
  "Macrophages"   = "#4DBBD5FF",
  "OL-mature"        = "#00A087FF",
  "Pericyte"       = "#3C5488FF",
  "OPCs"   = "#F39B7FFF",
  "other"   = "#8491B4FF"
)
ggplot(data, aes(x = 2, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", size = 0.6) +
  coord_polar(theta = "y", start = 0) +

  #geom_text(aes(x = 2.3,
  #              label = percent_label,
  #              angle = 90 - 180 * mid_angle / pi),
  #          size = 4.5,
  #          fontface = "bold",
  #          color = "black") +
  scale_fill_manual(
    name = "",
    values = colors[levels(data$celltype)],
    breaks = levels(data$celltype),
    labels = paste0(levels(data$celltype), " (",
                    sprintf("%.2f%%", data$percentage[match(levels(data$celltype), data$celltype)]), ")")
  ) +
  xlim(0.5, 2.5) +
  labs(
    title = "SPP1_POS Positive Cells",
    subtitle = paste0("Total cells: ", sum(data$count))
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40", margin = margin(b = 15)),
    legend.position = "right",
    legend.text = element_text(size = 10.5),
    legend.key.size = unit(0.7, "cm"),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(filename = "figs/SPP1.donut_chart.pdf",width = 8,height = 6)

mg=readRDS("../analysis_batch05/MG.cluster.v1.rds")
obj$subtype=NULL
mg$subtype=paste0("MG",as.numeric(mg$seurat_clusters))
colnames(mg@meta.data)
obj=AddMetaData(obj,mg@meta.data[,c("subtype"),drop=F])
p1 <- DimPlot(obj, reduction = "umap", group.by='subtype',label=T,label.size = 5)+NoLegend()
```

---

## 4. Neuron Sub-clustering

```r
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
library(tidyverse)

rawobj=readRDS("neuron.v1.rds")
p1 <- DimPlot(rawobj, reduction = "umap", group.by='celltype',label=T,label.size = 5)#+xlim(-15,20)
p2 <- DimPlot(rawobj, reduction = "umap", group.by='seurat_clusters',label=T)#+xlim(-15,20)
obj=subset(rawobj,celltype=="Non_neuronal",invert=T)
DefaultAssay(obj)<-'RNA'
obj <- SCTransform(obj, vst.flavor = "v2",variable.features.n=2000)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:20, reduction = "pca")
obj <- RunUMAP(obj, dims = 1:20, reduction = "pca",verbose =F)
obj <- FindClusters(obj, resolution = c(0.1,0.2,0.3,0.4,0.5),verbose = F)
p1<-DimPlot(obj, reduction = 'umap', group.by = "SCT_snn_res.0.1",label=T,label.size = 8)
p2<-DimPlot(obj, reduction = 'umap', group.by = "SCT_snn_res.0.2",label=T,label.size = 8)
p3<-DimPlot(obj, reduction = 'umap', group.by = "SCT_snn_res.0.3",label=T,label.size = 8)
p1<-DimPlot(obj, reduction = 'umap', group.by = "SCT_snn_res.0.4",label=T,label.size = 8)
p2<-DimPlot(obj, reduction = 'umap', group.by = "SCT_snn_res.0.5",label=T,label.size = 8)
markers=c(
"Slc17a7","Slc17a6","Camk2a","Neurod6","Tbr1","Cux2","Satb2","Rorb","Fezf2","Ctip2",
"Gad1",      # GAD67
"Gad2",      # GAD65
"Slc32a1",   # VGAT
"Pvalb","Sst","Vip","Lamp5","Npy","Cck","Cxcl14",
"Th","Slc6a3",    # DAT
"Foxa2","Nurr1",     # Nr4a2
"Tph2","Slc6a4",    # SERT
"Chat","Slc5a7","Th","Dbh","Slc6a2",    # NET
"Hdc","Tac1"
)
FeaturePlot(obj, features = unique(markers),ncol=7,label=F,pt.size=0.3,reduction='umap')
obj <- FindClusters(obj, resolution = 0.2,verbose = F)
p1<-DimPlot(obj, reduction = 'umap', group.by = "seurat_clusters",label=T,label.size = 8)
p2<-DimPlot(obj, reduction = 'umap', group.by = "Group",label=T,label.size = 4)
jjDotPlot(object = obj,
          gene = unique(markers),
          #xtree = T,
          dot.col = c('blue','white','red'),
          #rescale = T,
          rescale.min = -2,
          rescale.max = 2,
          midpoint = 0
         )+ theme(
    plot.margin = margin(r=0),
    axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5, face = "bold",size=12),
    axis.text.y = element_text(face = "bold",size=12),
    axis.title = element_blank(),
    legend.title = element_text(size = 12,face = "bold"),
    legend.text = element_text(size = 12) ,
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank()
  )#+coo

DefaultAssay(obj)="RNA"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,features=rownames(obj))
clustermarkers <- FindAllMarkers(object = obj, only.pos = T,verbose=F)
write.csv(clustermarkers,'MG.cluster.markers.csv')

cDEG=clustermarkers%>%
    filter(p_val_adj<0.001,pct.1>0.5)%>%
    arrange(cluster,desc(avg_log2FC))%>%
    filter(!str_detect(gene,'Gm|Rik|Rp'))|>
    #filter(gene%in%rownames(obj@assays$SCT$scale.data))|>
    group_by(cluster)|>
    top_n(n=20,wt=avg_log2FC)

p<-DoHeatmap(object = subset(obj, downsample = 1000),features = cDEG$gene)
split(cDEG$gene,cDEG$cluster)

library(scales)
fig=QCplot(obj,group.by="seurat_clusters")
fig1<-plot_cluster_number(obj,cluster_col="seurat_clusters",type="num")
fig2<-plot_cluster_number(obj,cluster_col="seurat_clusters",type="prop")
cluster_annotations <- c(
  "0" = "D2 MSN",
  "1" = "D1 MSN",
  "2" = "SST/PVALB interneuron",
  "3" = "GABAergic cortical neuron",
  "4" = "Glutamatergic neuron",
  "5" = "Deep layer neuron",
  "6" = "Glutamatergic neuron",
  "7" = "Immature neuron",
  "8" = "Glutamatergic neuron",
  "9" = "GABAergic cortical neuron",
  "10" = "Neural progenitor",
  "11" = "Npy/Sst interneuron",
  "12" = "Gata3+ GABAergic neuron"
)
obj$celltype <- unname(cluster_annotations[as.character(obj$seurat_clusters)])
obj$celltype <- factor(obj$celltype)
Idents(obj) <- obj$celltype
p1<-DimPlot(obj, group.by = "celltype", label = TRUE)
p2<-DimPlot(obj, group.by = "seurat_clusters", label = TRUE)
marker_genes <- list(
  "D2 MSN" = c("Drd2", "Adora2a", "Penk", "Gpr88"),
  "D1 MSN" = c("Drd1", "Pdyn", "Tac1", "Ebf1"),
  "SST/PVALB interneuron" = c("Sst", "Pvalb", "Avp", "Kcng4"),
  "GABAergic cortical neuron" = c("Gabrg1", "Aldoc", "Rmst", "Glra3"),
  "Glutamatergic neuron" = c("Slc17a6", "Lef1", "Zic5", "Lhx9"),
  "Deep layer neuron" = c("Erbb4", "Foxp2", "Pcdh8", "Chst9"),
  "Cck+ glutamatergic neuron" = c("Slc17a7", "Cck", "Rorb", "Prkcd"),
  "Immature neuron" = c("Neurod2", "Bhlhe22", "Nptx1"),
  "Stress-related neuron" = c("Atf3", "Hspb1", "Flnc", "Adcyap1"),
  "MGE/CGE interneuron" = c("Lhx6", "Maf", "Mafb", "Kit"),
  "Neural progenitor" = c("Sox4", "Sox11", "Dcx", "Tcf4"),
  "Npy/Sst interneuron" = c("Npy", "Nos1", "Sox6", "Sst"),
  "Gata3+ GABAergic neuron" = c("Gata3", "Tal1", "Chrna6", "Otx2")
)
all_markers_expanded <- unique(unlist(marker_genes))

#print(core_markers)
jjDotPlot(object = obj,
          gene = all_markers_expanded,
          #xtree = T,
          dot.col = c('blue','white','red'),
          #rescale = T,
          rescale.min = -2,
          rescale.max = 2,
          midpoint = 0
         )+ theme(
    plot.margin = margin(r=0),
    axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5, face = "bold",size=12),
    axis.text.y = element_text(face = "bold",size=12),
    axis.title = element_blank(),
    legend.title = element_text(size = 12,face = "bold"),
    legend.text = element_text(size = 12) ,
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank()
  )+NoLegend()
fig1<-plot_cluster_number(obj,cluster_col="celltype",type="num")
fig2<-plot_cluster_number(obj,cluster_col="celltype",type="prop")
excit_markers <- c("Slc17a6", "Slc17a7")
inhib_markers <- c("Gad1", "Gad2", "Slc32a1")
FeaturePlot(obj, features = c("Slc17a6", "Slc17a7", "Gad1", "Gad2"),
            ncol = 4, reduction = "umap")
DotPlot(obj, features = c(excit_markers, inhib_markers), group.by = "seurat_clusters") +
  RotatedAxis() +
  ggtitle("Excitatory vs Inhibitory Marker Validation")
excit_inhib_annotations <- c(
  "0" = "Inhibitory",
  "1" = "Inhibitory",
  "2" = "Inhibitory",
  "3" = "Inhibitory",
  "4" = "Excitatory",
  "5" = "Inhibitory",
  "6" = "Excitatory",
  "7" = "Excitatory",
  "8" = "Excitatory",
  "9" = "Inhibitory",
  "10" = "Inhibitory",
  "11" = "Inhibitory",
  "12" = "Inhibitory"
)
obj$celltype <- unname(excit_inhib_annotations[as.character(obj$seurat_clusters)])
obj$celltype <- factor(obj$celltype)
Idents(obj) <- obj$celltype
p1<-DimPlot(obj, group.by = "celltype", label = TRUE)
p2<-DimPlot(obj, group.by = "seurat_clusters", label = TRUE)
fig1<-plot_cluster_number(obj,cluster_col="celltype",type="num")
fig2<-plot_cluster_number(obj,cluster_col="celltype",type="prop")
saveRDS(obj,'neuron.v2.rds')
```

---

## 5. Neuron DEA (MAST)

```r
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)

obj=readRDS("neuron.v2.rds")
obj$Group=obj$Sample
obj$Group=ifelse(obj$Group=="CPSP-day7",'WT-TH',obj$Group)
obj$Group=ifelse(obj$Group=="saline",'WT-saline',obj$Group)
obj$Group=ifelse(obj$Group=="SPP1-ko-cpsp",'Spp1-ko-TH',obj$Group)
obj$Group=ifelse(obj$Group=="SPP1-ko-control",'Spp1-ko-saline',obj$Group)
DefaultAssay(obj)<-'RNA'
obj$main_type="neuron"

run_mast<-function(group1,group2){
    results <- FindMarkersAllCellTypes(
        seurat_obj = obj,
        group1 = group1,
        group2 = group2,
        group.by = "Group",
        celltype_col = "celltype",
        test.use="MAST",
        return_all = T,
        verbose=F)
    res1=results$combined_results
    results2 <- FindMarkersAllCellTypes(
        seurat_obj = obj,
        group1 = group1,
        group2 = group2,
        group.by = "Group",
        celltype_col = "main_type",
        test.use="MAST",
        return_all = T,
        verbose=F)
    res2=results2$combined_results
    res=rbind(res1,res2)
}

res1=run_mast("WT-saline","WT-TH")
res2=run_mast("Spp1-ko-saline","Spp1-ko-TH")
res3=run_mast("WT-TH","Spp1-ko-TH")
res=rbind(res1,res2,res3)
write.csv(res,'neuron.DEA.mast.raw.v2.csv',quote=F)
res=read.csv("neuron.DEA.mast.raw.v2.csv",row.names = 1)

ggVolcano<-function(result,padj.t=0.05,logFC.t=1,showgene=NULL,xMax=NA,yMax=NA){
    df=result%>%
        mutate(log10padj =-log10(result$padj))%>%
        mutate(FC_abs =abs(result$log2FoldChange))%>%
        arrange(desc(FC_abs))#%>%
        #drop_na()
    yintercept=log10(padj.t)*-1
    xintercept=logFC.t

     if ("Down regulated" %in% df$bias){
        cols=c("#D62728FF","#1F77B4FF","#999999")
    }else{
        cols=c("#D62728FF","#999999")
    }

    if (is.na(yMax)){
        min_padj=min(df$padj[!is.na(df$padj)])
        if (min_padj==0){
            yMax=100
        }else{
            yMax=max(df$log10padj)+5
        }
    }

    if (is.na(xMax)){
        xMax=max(abs(df$log2FoldChange))+1
    }
    counts <- table(df$bias)
    count_vec <- as.character(counts[as.character(df$bias)])
    df$bias <- paste0(df$bias, " (", count_vec, ")")
    vec=unique(df$bias)
    weights <- ifelse(grepl("Up", vec), 1,
                 ifelse(grepl("Down", vec), 2, 3))
    df$bias=factor(df$bias,levels=vec[order(weights)])
    title=df$Contrast[1]

    if (!is.null(showgene)){
        df=df|>
            mutate(highlight=case_when(
            gene_name%in%showgene ~ gene_name,
            TRUE ~ ""))
        }
    fig<-ggscatter(df,
              x = "log2FoldChange",
              y = "log10padj",
              ylim=c(0,yMax), xlim=c(-xMax,xMax),
              ylab = "-log10(padj)",
              title = title,
              color = "bias",label=NULL,
              size = 2,
              show.legend.text = F,
              palette =cols)+
              geom_hline(yintercept=yintercept,col='grey',linetype=5)+
              geom_vline(xintercept=c(xintercept*-1,xintercept),col='grey',linetype=5)+
    theme_bw()+
    theme(
    #legend.position =c(0.22,0.88),legend.background = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 12,face = "bold"),
        plot.title = element_text(hjust = 0.5,size=15,face = "bold"),
        axis.text.x=element_text(size=12,face="plain"),
        axis.text.y=element_text(size=12,face="plain"),
        axis.title.x=element_text(size=15,face="bold"),
        axis.title.y=element_text(size=15,face="bold"),
        aspect.ratio = 1
        )+
      labs(color=paste0('Padj <',padj.t,';|log2FC|>',logFC.t))+
        scale_color_manual(
        values = c("#D62728FF", "#1F77B4FF", "gray"),
        breaks = c(levels(df$bias)[0], levels(df$bias)[1],levels(df$bias)[2])
        )

    if (!is.null(showgene)){
        fig=fig+
                geom_text_repel(
            aes(label = highlight),
            size = 5,
            min.segment.length = 0.2,
            seed = 42,
            box.padding = 0.3,
            max.overlaps = Inf,
            #arrow = arrow(length = unit(0.01, "npc")),
            nudge_x = .1,
            nudge_y = .25,
            color = "black")
    }
}

plot_Volcano<-function(df,adjustP.cutoff=0.05,log2FC.cutoff=1,use.pvalue=T){
    if (use.pvalue) {
        df$padj=df$p_val
    }else{
        df$padj=df$p_val_adj
    }
    df$log2FoldChange=df$avg_log2FC
    df$bias="Not Changed"
    df$bias[df$padj < adjustP.cutoff & df$avg_log2FC >  log2FC.cutoff] = "Up regulated"
    df$bias[df$padj < adjustP.cutoff & df$avg_log2FC < -log2FC.cutoff] = "Down regulated"
    df$Contrast=df$comparison
    plot_list=list()
    for (type in unique(df$Contrast)){
        sub_df=df[df$Contrast==type,]
        plot_list[[type]]=ggVolcano(sub_df,padj.t = adjustP.cutoff, logFC.t = log2FC.cutoff)+ggtitle(type)+
    theme(
    plot.title = element_text(hjust = 0.5,vjust=0,size=20),
    legend.position = 'top'
    )+
  guides(color = guide_legend(title = NULL))
    }
    combined_plot <- wrap_plots(plot_list, ncol = 3, nrow = 1)
}
DEG=res|>#filter(p_val_adj<0.05)|>
    filter(celltype=="neuron")
plot_Volcano(DEG,adjustP.cutoff=0.05,log2FC.cutoff=0.58,use.pvalue=F)
res$entrez <-mapIds(x=org.Mm.eg.db,
                        keys = res$gene,
                        column = "ENTREZID",
                        keytype = "SYMBOL",
                        multiVals = "first")
res$signed_p <- -log10(res$p_val_adj) * sign(res$avg_log2FC)
sub_res=res|>filter(celltype=="neuron",p_val_adj<0.05,abs(avg_log2FC)>0.58)|>arrange(comparison)
sub_res$comparison=factor(sub_res$comparison,levels=c("WT-saline_vs_WT-TH",'Spp1-ko-saline_vs_Spp1-ko-TH','WT-TH_vs_Spp1-ko-TH'))
sub_res$bias=ifelse(sub_res$avg_log2FC>0,'up','down')
sub_res$comparison_bias=paste(sub_res$comparison,sub_res$bias,sep='.')
DEG_list=split(sub_res$gene,sub_res$comparison)
DEG_list2=split(sub_res$gene,sub_res$comparison_bias)
saveRDS(DEG_list,'results20260605/DEG_list.rds')

run_enrich<-function(gene_list,type="GO"){
      degs <- lapply(gene_list, FUN = function(x) {
        out <- suppressMessages(suppressWarnings(
          bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
        ))
     })

    if (type=="GO"){
        clusterplot <- compareCluster(geneCluster = degs, fun = "enrichGO", OrgDb = "org.Mm.eg.db",ont = "ALL", pvalueCutoff  = 1,qvalueCutoff  = 1)
        clusterplot <- setReadable(clusterplot, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
    }else{
        clusterplot <- compareCluster(geneCluster = degs, fun = "enrichKEGG",organism="mmu", pvalueCutoff=0.1)
        clusterplot <- setReadable(clusterplot, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
        clusterplot@compareClusterResult$Description <- gsub(pattern = " - Mus musculus (house mouse)", replacement = "", clusterplot@compareClusterResult$Description, fixed = T)
    }
}

split_compareCluster_to_enrichList <- function(xx) {
  df <- xx@compareClusterResult
  clusters <- unique(df$Cluster)
  clusters <- clusters[!is.na(clusters)]
  enrich_list <- list()

  for(cluster in clusters) {
    df_sub <- df[df$Cluster == cluster, ]

    if(nrow(df_sub) == 0) next
    genes_cluster <- character(0)
    if(!is.null(xx@geneClusters) && cluster %in% names(xx@geneClusters)) {
      genes_cluster <- as.character(xx@geneClusters[[cluster]])
    }
    enrich_obj <- new("enrichResult",
                      result = df_sub,
                      gene = genes_cluster,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      universe = character(0),
                      geneSets = list(),
                      organism = "Unknown",
                      keytype = if(!is.null(xx@keytype)) xx@keytype else "UNKNOWN",
                      ontology = "Unknown",
                      readable = if(!is.null(xx@readable)) xx@readable else FALSE)
    enrich_list[[cluster]] <- enrich_obj
  }
  return(enrich_list)
}

library(enrichplot)
library(DOSE)
GOres=run_enrich(DEG_list,type="GO")
KEGGres=run_enrich(DEG_list,type="KEGG")
GOdf=as.data.frame(GOres)
KEGGdf=as.data.frame(KEGGres)

#GOdf|>
#    filter(str_detect(Description,'ion tr'),p.adjust<0.01)
dotplot(GOres,showCategory=30)+  scale_y_discrete(labels=function(x) str_wrap(x, width=100))+ggtitle("GO")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
enrich_list=split_compareCluster_to_enrichList(GOres)
plots <- list()
for(name in names(enrich_list)) {
  p <- dotplot(enrich_list[[name]],
               showCategory = 20,
               title = name) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+  scale_y_discrete(labels=function(x) str_wrap(x, width=60))
  plots[[name]] <- p
}
combined <- ggarrange(plotlist = plots,
                      ncol = 3,
                      nrow = 1,
                      common.legend = TRUE,
                      legend = "left",
                      align = "v")
print(combined)
go_ids <- c(
  "GO:0042391",
  "GO:0010959",
  "GO:0034765",
  "GO:0048167",
  "GO:0006816",
  "GO:0060078",
  "GO:0055074",
  "GO:0001505",
  "GO:0097553",
  "GO:1900271",
  "GO:0034765",
  "GO:0010959",
  "GO:0042391",
  "GO:0001505",
  "GO:0048167",
  "GO:0060078",
  "GO:0006816",
  "GO:1900271",
  "GO:0055074"
)
GOres2=GOres|>
    filter(ID%in%go_ids,Cluster!="WT-TH_vs_Spp1-ko-TH")
GOdf2=as.data.frame(GOres2)

library(dplyr)
res_df <- as.data.frame(GOres2)
print(unique(res_df$Cluster))
res_df$neg_log10_p <- -log10(res_df$p.adjust)

if(!"Count" %in% colnames(res_df)) {
  if(is.character(res_df$GeneRatio)) {
    res_df$Count <- sapply(res_df$GeneRatio, function(x) {
      as.numeric(strsplit(x, "/")[[1]][1])
    })
  }
} else {
  res_df$Count <- as.numeric(res_df$Count)
}
res_df <- res_df %>%
  group_by(Cluster) %>%
  mutate(Description = reorder(Description, Count))
cluster_values <- unique(res_df$Cluster)
shape_map <- c(16, 15)
names(shape_map) <- cluster_values
p2 <- ggplot(res_df, aes(x = Count,
                          y = Description,
                          shape = Cluster,
                          color = neg_log10_p,
                          size = Count)) +
  geom_point(stroke = 1.2, alpha = 0.9) +
  scale_shape_manual(values = shape_map) +
  scale_color_gradient(low = "#327EBA", high = "#E06663",
                       name = expression(-log[10]("p.adjust"))) +
  scale_size_continuous(range = c(3, 10),
                        name = "Gene Count") +
 scale_x_continuous(
    limits = c(0, 50),
    expand = expansion(mult = c(0, 0.05))
  )+
 #xlim(1,50)+
  labs(x = "Gene Count",
       y = "GO Term",
       title = "GO Enrichment Comparison",
       shape = "Group") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray50"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())
ggsave(filename = "figs/neuron.DEG.GO.pdf",plot = p2,width = 10,height = 6)

GOres3=run_enrich(DEG_list2,type="GO")
GOdf3=as.data.frame(GOres3)
filename=get_outname('./',"neuron.up.down",'GO.tsv')
write_tsv(GOdf3,filename)
GOres4=GOres3|>
    filter(ID%in%go_ids)|>
    filter(!str_detect(Cluster,"WT-TH_vs_Spp1-ko-TH"))|>
    filter(str_detect(Cluster,"up"))
res_df <- as.data.frame(GOres4)
res_df$neg_log10_p <- -log10(res_df$p.adjust)

if(!"Count" %in% colnames(res_df)) {
  if(is.character(res_df$GeneRatio)) {
    res_df$Count <- sapply(res_df$GeneRatio, function(x) {
      as.numeric(strsplit(x, "/")[[1]][1])
    })
  }
} else {
  res_df$Count <- as.numeric(res_df$Count)
}
res_df <- res_df %>%
  group_by(Cluster) %>%
  mutate(Description = reorder(Description, Count))
cluster_values <- unique(res_df$Cluster)
shape_map <- c(15,16)
names(shape_map) <- cluster_values
p2 <- ggplot(res_df, aes(x = Count,
                          y = Description,
                          shape = Cluster,
                          color = neg_log10_p,
                          size = Count)) +
  geom_point(stroke = 1.2, alpha = 0.9) +
  scale_shape_manual(values = shape_map) +
  scale_color_gradient(low = "#327EBA", high = "#E06663",
                       name = expression(-log[10]("p.adjust"))) +
  scale_size_continuous(range = c(3, 10),
                        name = "Gene Count") +
  labs(x = "Gene Count",
       y = "GO Term",
       title = "GO Enrichment Comparison",
       shape = "Group") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray50"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())
print(p2)
GOres4=GOres3|>
    filter(ID%in%go_ids)|>
    filter(!str_detect(Cluster,"WT-TH_vs_Spp1-ko-TH"))|>
    filter(str_detect(Cluster,"down"))
res_df <- as.data.frame(GOres4)
print(unique(res_df$Cluster))
res_df$neg_log10_p <- -log10(res_df$p.adjust)

if(!"Count" %in% colnames(res_df)) {
  if(is.character(res_df$GeneRatio)) {
    res_df$Count <- sapply(res_df$GeneRatio, function(x) {
      as.numeric(strsplit(x, "/")[[1]][1])
    })
  }
} else {
  res_df$Count <- as.numeric(res_df$Count)
}
res_df <- res_df %>%
  group_by(Cluster) %>%
  mutate(Description = reorder(Description, Count))
cluster_values <- unique(res_df$Cluster)
shape_map <- c(15,16)
names(shape_map) <- cluster_values
p2 <- ggplot(res_df, aes(x = Count,
                          y = Description,
                          shape = Cluster,
                          color = neg_log10_p,
                          size = Count)) +
  geom_point(stroke = 1.2, alpha = 0.9) +
  scale_shape_manual(values = shape_map) +
  scale_color_gradient(low = "#327EBA", high = "#E06663",
                       name = expression(-log[10]("p.adjust"))) +
  scale_size_continuous(range = c(3, 10),
                        name = "Gene Count") +
  labs(x = "Gene Count",
       y = "GO Term",
       title = "GO Enrichment Comparison",
       shape = "Group") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray50"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())
print(p2)
dotplot(KEGGres,showCategory=30)+  scale_y_discrete(labels=function(x) str_wrap(x, width=100))+ggtitle("GO")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
enrich_list=split_compareCluster_to_enrichList(KEGGres)
plots <- list()
for(name in names(enrich_list)) {
  p <- dotplot(enrich_list[[name]],
               showCategory = 20,
               title = name) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+  scale_y_discrete(labels=function(x) str_wrap(x, width=60))
  plots[[name]] <- p
}

outdir="./results20260605/"
filename=get_outname(outdir,"neuron",'GO.tsv')
write_tsv(GOdf,filename)
filename=get_outname(outdir,'neuron','KEGG.tsv')
write_tsv(KEGGdf,filename)
```

---

## 6. Neuron DEA (Wilcoxon)

```r
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)

obj=readRDS("neuron.v2.rds")
obj$Group=obj$Sample
obj$Group=ifelse(obj$Group=="CPSP-day7",'WT-TH',obj$Group)
obj$Group=ifelse(obj$Group=="saline",'WT-saline',obj$Group)
obj$Group=ifelse(obj$Group=="SPP1-ko-cpsp",'Spp1-ko-TH',obj$Group)
obj$Group=ifelse(obj$Group=="SPP1-ko-control",'Spp1-ko-saline',obj$Group)
DefaultAssay(obj)<-'RNA'
obj$main_type="neuron"

run_mast<-function(group1,group2){
    results <- FindMarkersAllCellTypes(
        seurat_obj = obj,
        group1 = group1,
        group2 = group2,
        group.by = "Group",
        celltype_col = "celltype",
        test.use="wilcox",
        return_all = T,
        verbose=F)
    res1=results$combined_results
    results2 <- FindMarkersAllCellTypes(
        seurat_obj = obj,
        group1 = group1,
        group2 = group2,
        group.by = "Group",
        celltype_col = "main_type",
        test.use="wilcox",
        return_all = T,
        verbose=F)
    res2=results2$combined_results
    res=rbind(res1,res2)
}

res1=run_mast("WT-saline","WT-TH")
res2=run_mast("Spp1-ko-saline","Spp1-ko-TH")
res3=run_mast("WT-TH","Spp1-ko-TH")
res=rbind(res1,res2,res3)
write.csv(res,'neuron.DEA.wilcox.raw.v2.csv',quote=F)

plot_jjVolcano<-function(df,adjustP.cutoff=0.05,log2FC.cutoff=1,topn=5){
    df$cluster=df$celltype
    df$avg_log2FC=ifelse(df$avg_log2FC>10,10,df$avg_log2FC)
    df$avg_log2FC=ifelse(df$avg_log2FC< -10,-10,df$avg_log2FC)
    Contrast=df$comparison
    cluster_count=length(unique(df$cluster))
    p<-jjVolcano(diffData = df,
          log2FC.cutoff=log2FC.cutoff,
          adjustP.cutoff=adjustP.cutoff,
          topGeneN = topn,
          tile.col=paletteer::paletteer_c("grDevices::Dark 2", cluster_count),
          base_size=16,
         )+ggtitle(Contrast)+
    theme(
    plot.title = element_text(hjust = 0.5,vjust=0,size=20),
    legend.position = 'none'
    )
}

plot_Volcano<-function(df,adjustP.cutoff=0.05,log2FC.cutoff=1,use.pvalue=T){
    if (use.pvalue) {
        df$padj=df$p_val
    }else{
        df$padj=df$p_val_adj
    }
    df$log2FoldChange=df$avg_log2FC
    df$bias="Not Changed"
    df$bias[df$padj < adjustP.cutoff & df$avg_log2FC >  log2FC.cutoff] = "Up regulated"
    df$bias[df$padj < adjustP.cutoff & df$avg_log2FC < -log2FC.cutoff] = "Down regulated"
    df$Contrast=paste0(df$comparison,':',df$celltype)
    plot_list=list()
    for (type in unique(df$Contrast)){
        sub_df=df[df$Contrast==type,]
        plot_list[[type]]=ggVolcano(sub_df,padj.t = adjustP.cutoff, logFC.t = log2FC.cutoff)+ggtitle(type)
    }
    combined_plot <- wrap_plots(plot_list, ncol = 3, nrow = 1)
}
DEG=res|>#filter(p_val_adj<0.05)|>
    filter(comparison=="WT-saline_vs_WT-TH")
plot_jjVolcano(DEG,adjustP.cutoff=0.05,log2FC.cutoff=0.58,topn=5)
DEG=res|>#filter(p_val_adj<0.05)|>
    filter(comparison=="WT-saline_vs_WT-TH")
plot_Volcano(DEG)
DEG=res|>#filter(p_val_adj<0.01)|>
    filter(comparison=="Spp1-ko-saline_vs_Spp1-ko-TH")
plot_jjVolcano(DEG,adjustP.cutoff=0.05,log2FC.cutoff=0.58,topn=5)
DEG=res|>#filter(p_val_adj<0.05)|>
    filter(comparison=="Spp1-ko-saline_vs_Spp1-ko-TH")
plot_Volcano(DEG)
DEG=res|>#filter(p_val_adj<0.01)|>
    filter(comparison=="WT-TH_vs_Spp1-ko-TH")
plot_jjVolcano(DEG,adjustP.cutoff=0.05,log2FC.cutoff=0.58,topn=5)
DEG=res|>#filter(p_val_adj<0.05)|>
    filter(comparison=="WT-TH_vs_Spp1-ko-TH")
plot_Volcano(DEG)
res$entrez <-mapIds(x=org.Mm.eg.db,
                        keys = res$gene,
                        column = "ENTREZID",
                        keytype = "SYMBOL",
                        multiVals = "first")
res$signed_p <- -log10(res$p_val_adj) * sign(res$avg_log2FC)

contrasts <- c("WT-saline_vs_WT-TH",
             "Spp1-ko-saline_vs_Spp1-ko-TH",
             "WT-TH_vs_Spp1-ko-TH")
outdir <- "results20260605"
for (Contrast in contrasts) {
  sub_res <- res |> filter(comparison == Contrast)
  DEG_list <- split(sub_res, sub_res$celltype)
  DEG_list <- lapply(DEG_list, function(x) {
    x <- x |>
      filter(!is.na(entrez)) |>
      distinct(entrez, .keep_all = T) |>
      arrange(desc(avg_log2FC))
    fc_vector <- x$avg_log2FC
    names(fc_vector) <- as.character(x$entrez)
  })
  GOclusterplot <- compareCluster(geneCluster = DEG_list, fun = "gseGO",
                                  OrgDb = "org.Mm.eg.db", ont = "ALL")
  GOclusterplot <- setReadable(GOclusterplot, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  GOdf <- as.data.frame(GOclusterplot)
  print(dotplot(GOclusterplot, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign) +
          scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
          theme(axis.text.x = element_text(angle = 60, hjust = 1, colour = "black", size = 12)) +
          ggtitle(paste("GO -", Contrast)))
  KEGGclusterplot <- compareCluster(geneCluster = DEG_list, fun = "gseKEGG",
                                    organism = "mmu", pvalueCutoff = 0.05)
  KEGGclusterplot <- setReadable(KEGGclusterplot, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  KEGGclusterplot@compareClusterResult$Description <-
    gsub(" - Mus musculus (house mouse)", "", KEGGclusterplot@compareClusterResult$Description, fixed = T)
  KEGGdf <- as.data.frame(KEGGclusterplot)
  print(dotplot(KEGGclusterplot, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign) +
          scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
          theme(axis.text.x = element_text(angle = 60, hjust = 1, colour = "black", size = 12)) +
          ggtitle(paste("KEGG -", Contrast)))

  write_tsv(GOdf,   get_outname(outdir, Contrast, "GO.tsv"))
  write_tsv(KEGGdf, get_outname(outdir, Contrast, "KEGG.tsv"))
}
```

## 7. MG Clustering

```r
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
library(tidyverse)
rawobj=readRDS("scRNA.celltype.rds")
p1 <- DimPlot(rawobj, reduction = "umap", group.by='celltype',label=T,label.size = 5)#+xlim(-15,20)
p2 <- DimPlot(rawobj, reduction = "umap", group.by='seurat_clusters',label=T)#+xlim(-15,20)
obj=subset(rawobj,celltype=="Microglia")
DefaultAssay(obj)<-'RNA'
obj <- SCTransform(obj, vst.flavor = "v2",variable.features.n=2000)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:20, reduction = "pca")
obj <- RunUMAP(obj, dims = 1:20, reduction = "pca",verbose =F)
obj <- FindClusters(obj, resolution = c(0.1,0.2,0.3,0.4,0.5),verbose = F)
p1<-DimPlot(obj, reduction = 'umap', group.by = "SCT_snn_res.0.1",label=T,label.size = 8)
p2<-DimPlot(obj, reduction = 'umap', group.by = "SCT_snn_res.0.2",label=T,label.size = 8)
p3<-DimPlot(obj, reduction = 'umap', group.by = "SCT_snn_res.0.3",label=T,label.size = 8)
p1<-DimPlot(obj, reduction = 'umap', group.by = "SCT_snn_res.0.4",label=T,label.size = 8)
p2<-DimPlot(obj, reduction = 'umap', group.by = "SCT_snn_res.0.5",label=T,label.size = 8)
obj <- FindClusters(obj, resolution = 0.1)

obj=readRDS("MG.cluster.v1.rds")
obj$subtype=paste0("MG",as.numeric(obj$seurat_clusters))
p1 <- DimPlot(obj, reduction = "umap", group.by='subtype',label=T,label.size = 8)+NoLegend()
ggsave("figs/MG.cluster.pdf",p1)

obj$Group=obj$Sample
obj$Group=ifelse(obj$Group=="saline",'Saline',obj$Group)
obj$Group=ifelse(obj$Group=="CPSP-day7",'Th-day7',obj$Group)
obj$Group=ifelse(obj$Group=="CPSP-3w",'Th-day21',obj$Group)
obj$Group=ifelse(obj$Group=="Cx3crl-cOpn5-No-Tamox",'Opto-NoTAM',obj$Group)
obj$Group=ifelse(obj$Group=="cOpn5-day1",'Opto-day1',obj$Group)
obj$Group=ifelse(obj$Group=="Cx3crl-cOpn5-3w",'Opto-day21',obj$Group)
obj$Group=factor(obj$Group,levels=c("Saline",'Th-day7','Th-day21','Opto-NoTAM','Opto-day1','Opto-day21','SPP1-ko-control','SPP1-ko-cpsp'))
p1<-DimPlot(obj, reduction = 'umap', group.by = "subtype",label=F,label.size = 4,split.by="Group",ncol=3)
ggsave("figs/MG.sample.pdf",p1,width = 8,height = 8)

genelist=c('Cx3cr1','P2ry12','Tmem119',"Hexb","Siglech",'Iba1','Spp1')
FeaturePlot(obj, features = genelist ,label=T,ncol=3,reduction = 'umap')
ggplot(obj@meta.data, aes(x=Group, fill=Idents(obj))) +
        geom_bar(position = "fill",width = 0.7,size = 0.2,colour = '#222222')+
        theme_classic()+
        theme(
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5,size = 20),
        axis.text.y=element_text(size = 19),
        axis.text = element_text(color = 'black', size = 12),
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
             )+
        #scale_fill_manual(values = my.colors)+
        #scale_fill_gradient(low="red",high="blue")+
        labs(x=NULL,y="Percentage")
ggsave("figs/MG.subtype.percent.pdf",width = 6,height = 6)

test=subset(obj,Group%in%c("Saline",'Th-day7','Th-day21'))
FeaturePlot(test, features = "Spp1" ,label=T,ncol=3,reduction = 'umap',split.by='Group')+ theme(legend.position = "right")
ggsave("figs/MG.spp1.1.pdf",width = 15,height = 5)

test=subset(obj,Group%in%c('Opto-NoTAM','Opto-day1','Opto-day21'))
FeaturePlot(test, features = "Spp1" ,label=T,ncol=3,reduction = 'umap',split.by='Group')+ theme(legend.position = "right")
ggsave("figs/MG.spp1.2.pdf",width = 15,height = 5)

test=subset(obj,Group%in%c('SPP1-ko-control','SPP1-ko-cpsp'))
FeaturePlot(test, features = "Spp1" ,label=T,ncol=3,reduction = 'umap',split.by='Group')+ theme(legend.position = "right")
ggsave("figs/MG.spp1.3.pdf",width = 10,height = 5)

saveRDS(obj,'MG.cluster.v1.rds')

obj=readRDS("MG.cluster.v1.rds")
#obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = TRUE)

FeaturePlot(obj,features=c("Cxcl14",'Spp1','Cxcl10','Igf1'))
DefaultAssay(obj)<-"RNA"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,features=rownames(obj))
clustermarkers <- FindAllMarkers(object = obj, only.pos = T)
write.csv(clustermarkers,'MG.cluster.markers.csv')

cDEG=clustermarkers%>%
    filter(p_val_adj<0.001,pct.1>0.5)%>%
    arrange(cluster,desc(avg_log2FC))%>%
    filter(!str_detect(gene,'Gm|Rik|Rp'))|>
    #filter(gene%in%rownames(obj@assays$SCT$scale.data))|>
    group_by(cluster)|>
    top_n(n=20,wt=avg_log2FC)
clustermarkers|>
    filter(gene%in%c("P2ry12",'Siglech','Gpr34','Socs3','Hexb','Olfml3','Fcrls'))
p<-DoHeatmap(object = subset(obj, downsample = 1000),features = c(cDEG$gene,"P2ry12",'Siglech','Gpr34','Socs3','Hexb','Olfml3','Fcrls'))
p<-DoHeatmap(object = subset(obj, downsample = 1000),features = cDEG$gene)
ggsave("figs/MG.cluster.markers.pdf",width = 12,height = 12)

library(clusterProfiler)
library(org.Mm.eg.db)
cDEG2=clustermarkers%>%
    filter(p_val_adj<0.001,avg_log2FC>1,pct.1>0.5)%>%
    arrange(cluster,desc(avg_log2FC))%>%
    distinct(gene,.keep_all = T)
cDEG2$entrez <-mapIds(x=org.Mm.eg.db,
                        keys = cDEG2$gene,
                        column = "ENTREZID",
                        keytype = "SYMBOL",
                        multiVals = "first")
cDEG2<- subset(cDEG2, !is.na(entrez))
#cDEG2$cluster=paste0('cluster',as.character(cDEG2$cluster))
DEG_list <- split(cDEG2$entrez, cDEG2$cluster)
GOclusterplot <- compareCluster(geneCluster = DEG_list, fun = "enrichGO", OrgDb = "org.Mm.eg.db")
GOclusterplot <- setReadable(GOclusterplot, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
fig<-dotplot(GOclusterplot,showCategory = 10)+
    scale_y_discrete(labels=function(x) str_wrap(x, width=100))+
    theme(axis.text.x=element_text(angle=0,hjust = 0.5,colour="black",size=12))
ggsave("figs/MG.cluster.DEG.GO.pdf",width = 10,height = 8)
KEGGclusterplot <- compareCluster(geneCluster = DEG_list, fun = "enrichKEGG",organism="mmu", pvalueCutoff=0.05)
KEGGclusterplot <- setReadable(KEGGclusterplot, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
KEGGclusterplot@compareClusterResult$Description <- gsub(pattern = " - Mus musculus (house mouse)", replacement = "", KEGGclusterplot@compareClusterResult$Description, fixed = T)
dotplot(KEGGclusterplot,showCategory = 10)+
    scale_y_discrete(labels=function(x) str_wrap(x, width=100))+
    theme(axis.text.x=element_text(angle=0,hjust = 0.5,colour="black",size=12))
ggsave("figs/MG.cluster.DEG.KEGG.pdf",width = 9.5,height = 8)
```

## 8. MG CPSP DEG

```r
library(Seurat )
library(tidyverse)
library(gtools)
library(cowplot)
library(gridExtra)
library(patchwork)
library(scales)
library(viridis)
library(paletteer)
library(RColorBrewer)
library(ggpubr)
library(IRdisplay)
library(clusterProfiler)
library(org.Mm.eg.db)

library(ComplexHeatmap)
library(RColorBrewer)
options(warn=-1)
scRNA=readRDS("scRNA.celltype.rds")
scRNA@version
scRNA@meta.data$Sample=scRNA$Group
scRNA$Group=ifelse(str_detect(scRNA$Sample,"CPSP"),'CPSP',scRNA$Sample)
scRNA$Group=ifelse(str_detect(scRNA$Sample,"cOpn5"),'cOpn5',scRNA$Group)
scRNA$Group=ifelse(str_detect(scRNA$Sample,"SPP1"),'SPP1-KO',scRNA$Group)
scRNA$Group=ifelse(str_detect(scRNA$Sample,"ITGB1"),'ITGB1',scRNA$Group)
p1 <- DimPlot(scRNA, reduction = "umap", group.by='celltype',label=T,label.size = 5)
p2 <- DimPlot(microglia, reduction = "umap", group.by='seurat_clusters',label=T,label.size=8)

GO_info=read.csv("~/DATA/Ref/musculus/mm10.anno.tsv",sep='\t')
out=GO_info|>
    filter(str_detect(GOnames,'phagocytosis|engulfment|cytokine'))
GO_genelist=unique(out$gene_name)
length(GO_genelist)
GO_info2=read.csv("~/DATA/Ref/musculus/GO.csv")
out=GO_info2|>
    filter(str_detect(name,'phagocytosis|ECM'))
GO_ids=unique(out$id)
length(GO_ids)
p1 <- DimPlot(scRNA, reduction = "umap", group.by='celltype',label=F,label.size = 5,split.by="Sample",ncol=4)

get_DEG<-function(s.in,Sample1,Sample2,celltypelist){
    Contrast=paste(Sample1,Sample2,sep="_vs_")
    sub.in=subset(s.in,Sample%in%c(Sample1,Sample2))
    rawdegs.list = lapply(celltypelist, function(x){
        s.obj=subset(sub.in,celltype==x)
        out=FindMarkers(s.obj,ident.1 = Sample1, ident.2 = Sample2,logfc.threshold=0,min.pct = 0, recorrect_umi = FALSE)
        out$celltype=x
        out$gene=rownames(out)
    })
    sample.markers=do.call(rbind,lapply(rawdegs.list, function(x){
    x#table(x$avg_log2FC > 0 )
    }))
    sample.markers$cluster=sample.markers$celltype
    sample.markers$Contrast=Contrast
}
typelist=unique(scRNA$celltype)
typelist=typelist[!typelist%in%c("DCs",'Mast cells','T cells')]
res1=get_DEG(scRNA,"CPSP-day7",'saline',typelist)
res2=get_DEG(scRNA,"CPSP-3w",'saline',typelist)
res3=get_DEG(scRNA,"cOpn5-day1",'Cx3crl-cOpn5-No-Tamox',typelist)
res4=get_DEG(scRNA,'Cx3crl-cOpn5-3w',"Cx3crl-cOpn5-No-Tamox",typelist)
res5=get_DEG(scRNA,'SPP1-ko-cpsp',"SPP1-ko-control",typelist)
res=rbind(res1,res2,res3,res4,res5)

write.csv(res,'Sample.DEG.csv',quote=F)

plot_volcano<-function(data,cell,padj_threshold=0.001,logFC_threshold=1){
    df=data[data$celltype==cell,]
    label="gene"
    genes=c("Spp1","Fn1","Cxcl2","Ifi27l2a","Vim","Il1b","Apoe","Apoc1","Igf1","Cst7")
    genes=c(genes,GO_genelist)
    df$log10padj=-log10(df$p_val_adj)
    df$bias="Not Changed"
    df$bias[df$p_val_adj < padj_threshold & df$avg_log2FC >  logFC_threshold] = "Up regulated"
    df$bias[df$p_val_adj < padj_threshold & df$avg_log2FC < -logFC_threshold] = "Down regulated"
    df$bias=factor(df$bias,levels=c("Down regulated","Not Changed","Up regulated"))
    title=df$celltype[1]#paste(df$Contrast[1],df$celltype[1],sep=' , ')
    yintercept=log10(padj_threshold)*-1
    xintercept=logFC_threshold
    yMax=max(df$log10padj)+5
    xMax=max(df$avg_log2FC)+1
    df2=df[df$bias!="Not Changed",]
    showgene=genes[genes%in%df2$gene]

    if ("Down regulated" %in% df$bias){
        cols=c("#00AFBB","#999999","#FC4E07")
    }else{
        cols=c("#999999","#FC4E07")
    }
    fig<-ggscatter(df,
              x = "avg_log2FC",
              y = "log10padj",
              ylim=c(0,yMax), xlim=c(-xMax,xMax),
              ylab = "-log10(padj)",
              title = title,
              color = "bias",
              size = 1,
              label = label,
              #repel = T,
              show.legend.text = F,
              palette =cols ,
              label.select = showgene)+
              geom_hline(yintercept=yintercept,col='grey',linetype=5)+
              geom_vline(xintercept=c(xintercept*-1,xintercept),col='grey',linetype=5)+
              theme_bw()+theme(plot.title = element_text(hjust = 0.5))
    }

for (contrast in unique(res$Contrast)){
        data=res|>
            filter(Contrast==contrast)
        plot_list=lapply(typelist,function(cell) plot_volcano(data,cell))
        fig=wrap_plots(plot_list,guides = 'collect',ncol=9)+
            plot_annotation(contrast,theme=theme(plot.title=element_text(hjust=0,size=20)))#+
            #theme(legend.position = "top")
}
contrast="CPSP-day7_vs_saline"
data=res|>
        filter(Contrast==contrast)
p<-plot_volcano(data,'Microglia')+NoLegend()+ggtitle(contrast)
ggsave(paste('figs/',contrast,'.volcano.pdf'))

contrast="CPSP-day7_vs_saline"
data=res|>
        filter(Contrast==contrast)
p<-plot_volcano(data,'Microglia')+NoLegend()+ggtitle(contrast)
ggsave(paste('figs/',contrast,'.volcano.pdf'))

contrast="CPSP-3w_vs_saline"
data=res|>
        filter(Contrast==contrast)
p<-plot_volcano(data,'Microglia')+NoLegend()+ggtitle(contrast)
ggsave(paste('figs/',contrast,'.volcano.pdf'))

contrast="cOpn5-day1_vs_Cx3crl-cOpn5-No-Tamox"
data=res|>
        filter(Contrast==contrast)
p<-plot_volcano(data,'Microglia')+NoLegend()+ggtitle(contrast)
ggsave(paste('figs/',contrast,'.volcano.pdf'))

contrast="Cx3crl-cOpn5-3w_vs_Cx3crl-cOpn5-No-Tamox"
data=res|>
        filter(Contrast==contrast)
p<-plot_volcano(data,'Microglia')+NoLegend()+ggtitle(contrast)
ggsave(paste('figs/',contrast,'.volcano.pdf'))

res=read.csv("Sample.DEG.csv",row.names = 1)
DEG=res|>
    filter(cluster=="Microglia",p_val_adj<0.001,abs(avg_log2FC)>1)
DEG$entrez <-mapIds(x=org.Mm.eg.db,
                        keys = DEG$gene,
                        column = "ENTREZID",
                        keytype = "SYMBOL",
                        multiVals = "first")
DEG<- subset(DEG, !is.na(entrez))
DEG_list <- split(DEG$entrez, DEG$Contrast)
GOclusterplot <- compareCluster(geneCluster = DEG_list, fun = "enrichGO", OrgDb = "org.Mm.eg.db")
GOclusterplot <- setReadable(GOclusterplot, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
dotplot(GOclusterplot,showCategory = 10)+
    scale_y_discrete(labels=function(x) str_wrap(x, width=100))+
    theme(axis.text.x=element_text(angle=90,hjust = 1,colour="black",size=9))

ggsave("MG.DEG.enrich/microglia.GO.pdf",width = 10,height = 10)
Dir='MG.DEG.enrich'
GOdf=as.data.frame(GOclusterplot)
write_tsv(GOdf,paste0(Dir,'/',"MG",'.sampleDEG.GO.tsv'))
GOdf|>
    filter(str_detect(Description,"phagocytosis|ECM|Inter|interferon"))
KEGGclusterplot <- compareCluster(geneCluster = DEG_list, fun = "enrichKEGG",organism="mmu", pvalueCutoff=0.05)
KEGGclusterplot <- setReadable(KEGGclusterplot, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
KEGGclusterplot@compareClusterResult$Description <- gsub(pattern = " - Mus musculus (house mouse)", replacement = "", KEGGclusterplot@compareClusterResult$Description, fixed = T)
dotplot(KEGGclusterplot,showCategory = 10)+
    scale_y_discrete(labels=function(x) str_wrap(x, width=100))+
    theme(axis.text.x=element_text(angle=90,hjust = 1,colour="black",size=12))

ggsave("MG.DEG.enrich/microglia.KEGG.pdf",width = 10,height = 10)
KEGGdf=as.data.frame(KEGGclusterplot)
write_tsv(KEGGdf,paste0(Dir,'/',"MG",'.sampleDEG.KEGG.tsv'))
KEGGdf|>
    filter(str_detect(Description,"phagocytosis|ECM|Interferon|interferon"))
data=read.csv("Sample.DEG.csv",row.names = 1)
data$entrez <-mapIds(x=org.Mm.eg.db,
                        keys = data$gene,
                        column = "ENTREZID",
                        keytype = "SYMBOL",
                        multiVals = "first")
data<- subset(data, !is.na(entrez))|>arrange(desc(avg_log2FC))
sub_data=data|>filter(celltype=="Microglia")
DEG_list <- split(sub_data, sub_data$Contrast)
DEG_list <- lapply(DEG_list, function(x) {
  fc_vector <- x$avg_log2FC
  return(fc_vector)
})
GOclusterplot <- compareCluster(geneCluster = DEG_list, fun = "gseGO", OrgDb = "org.Mm.eg.db")
GOclusterplot <- setReadable(GOclusterplot, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
GOdf=as.data.frame(GOclusterplot)
dotplot(GOclusterplot, showCategory=10, split=".sign") + facet_grid(.~.sign)+
    scale_y_discrete(labels=function(x) str_wrap(x, width=100))+
    theme(axis.text.x=element_text(angle=60,hjust = 1,colour="black",size=12))
ggsave("MG.DEG.enrich/microglia.gsea.GO.pdf",width = 16,height = 20)

Dir='MG.DEG.enrich'
GOdf=as.data.frame(GOclusterplot)
write_tsv(GOdf,paste0(Dir,'/microglia.gsea.GO.tsv'))
KEGGclusterplot <- compareCluster(geneCluster = DEG_list, fun = "gseKEGG",organism="mmu", pvalueCutoff=0.05)
KEGGclusterplot <- setReadable(KEGGclusterplot, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
KEGGclusterplot@compareClusterResult$Description <- gsub(pattern = " - Mus musculus (house mouse)", replacement = "", KEGGclusterplot@compareClusterResult$Description, fixed = T)
KEGGdf=as.data.frame(KEGGclusterplot)

dotplot(KEGGclusterplot, showCategory=20, split=".sign") + facet_grid(.~.sign)+
    scale_y_discrete(labels=function(x) str_wrap(x, width=100))+
    theme(axis.text.x=element_text(angle=60,hjust = 1,colour="black",size=12))
ggsave("MG.DEG.enrich/microglia.gsea.KEGG.pdf",width = 16,height = 15)

KEGGdf=as.data.frame(KEGGclusterplot)
write_tsv(KEGGdf,paste0(Dir,'/microglia.gsea.KEGG.tsv'))
```

## 9. CellChat 

```r
library(ComplexHeatmap)
library(patchwork)
library(Seurat )
library(CellChat)
library(tidyverse)
library(IRdisplay)
scRNA=readRDS("scRNA.celltype.rds")
p1 <- DimPlot(scRNA, reduction = "umap", group.by='celltype',label=T,label.size = 5)#+xlim(-15,20)
p2 <- DimPlot(scRNA, reduction = "umap", group.by='Group',label=F)#+xlim(-15,20)
p1 + p2
subRNA=subset(scRNA,celltype%in%c("Astrocyte",'ECs','Macrophages','Microglia','Neuron','OL-diff','OL-mature','OPCs','Pericyte'))
DefaultAssay(subRNA)='RNA'
subRNA <- NormalizeData(subRNA)
CellChatDB <- CellChatDB.mouse
rm("scRNA")
s.list <- SplitObject(subRNA, split.by = "Group")
s.list <- lapply(X = s.list, FUN = function(x) {
    data.input = GetAssayData(object = x, assay = "RNA", slot = "data")
    meta.data =  x@meta.data
    meta.data = meta.data[!is.na(meta.data$celltype),]
    data.input = data.input[,row.names(meta.data)]
    meta.data$seurat_annotations = factor(meta.data$celltype)
    cellchat <- createCellChat(object = data.input,
                           meta = meta.data,
                           group.by = "celltype")
    CellChatDB.use <- CellChatDB # simply use the default CellChatDB
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat,only.pos = FALSE)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    #Compute the communication probability and infer cellular communication network
    cellchat <- computeCommunProb(cellchat,population.size = F)
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    #Infer the cell-cell communication at a signaling pathway level
    cellchat <- computeCommunProbPathway(cellchat)
    #Calculate the aggregated cell-cell communication network
    cellchat <- aggregateNet(cellchat)
})
saveRDS(s.list,'s.list.rds')


s.list=readRDS("s.list.rds")
object.list <- list(spp1_ko_cpsp = s.list[["SPP1-ko-cpsp"]], cpsp_day7 = s.list[["CPSP-day7"]])

for (i in 1:length(object.list)) {
    object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]],slot.name = "netP")
}
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
par(mfrow = c(1,2), xpd=TRUE)
gg1<-netVisual_diffInteraction(cellchat, weight.scale = T)
gg2<-netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
gg1 <- netVisual_heatmap(cellchat,font.size = 12,font.size.title = 15)
gg2 <- netVisual_heatmap(cellchat, measure = "weight",font.size = 12,font.size.title = 15)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(gg1, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(gg2, newpage = FALSE)
upViewport()
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

for (type in unique(subRNA$celltype)){
    gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = type,comparison = c(1, 2),label.size = 6,font.size.title = 20,font.size = 15)
    display(gg1)
}
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE, font.size = 12,thresh = 0.01)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE, font.size = 12,thresh = 0.01)
gg1 + gg2

# combining all the identified signaling pathways from different datasets
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i],
                                        width = 9, height = 20,font.size = 10, font.size.title = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1],
                                       width = 9, height = 20,font.size = 10, font.size.title = 15)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht1, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(ht2, newpage = FALSE)
upViewport()
gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:9),  comparison = c(1, 2), max.dataset = 2,font.size = 14,font.size.title = 12,
                        title.name = "Increased signaling in CPSP-day7", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:9),  comparison = c(1, 2), max.dataset = 1,font.size = 14,font.size.title = 12,
                        title.name = "Decreased signaling in CPSP-day7", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
df.net <- subsetCommunication(cellchat)
pathways.show <- c("SPP1")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
ht1<-netVisual_heatmap(object.list[[1]], signaling = pathways.show, color.heatmap = "Reds")+ggtitle(names(object.list[1]))
ht2<-netVisual_heatmap(object.list[[2]], signaling = pathways.show, color.heatmap = "Reds")+ggtitle(names(object.list[2]))
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht1, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(ht2, newpage = FALSE)
upViewport()
p1<-netAnalysis_contribution(object.list[[1]], signaling = pathways.show)+ggtitle(names(object.list[1]))
p2<-netAnalysis_contribution(object.list[[2]], signaling = pathways.show)+ggtitle(names(object.list[2]))
p1|p2
p1<-netVisual_bubble(object.list[[1]], sources.use = 4, targets.use = c(1:9), remove.isolate = FALSE,signaling =pathways.show,font.size = 20)+ggtitle(names(object.list[1]))
p2<-netVisual_bubble(object.list[[2]], sources.use = 4, targets.use = c(1:9), remove.isolate = FALSE,signaling =pathways.show,font.size = 20)+ggtitle(names(object.list[2]))
p1|p2
fig<-plotGeneExpression(cellchat, signaling = pathways.show, split.by = "datasets", colors.ggplot = T, type = "violin")
pathways.show <- c("CCL")

pdf("figs/cpsp_day7_vs_spp1_ko_cpsp.cellchat.1.pdf",width=12,height = 6)
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

dev.off()
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
ht1<-netVisual_heatmap(object.list[[1]], signaling = pathways.show, color.heatmap = "Reds")+ggtitle(names(object.list[1]))
ht2<-netVisual_heatmap(object.list[[2]], signaling = pathways.show, color.heatmap = "Reds")+ggtitle(names(object.list[2]))
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht1, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(ht2, newpage = FALSE)
upViewport()

pdf("figs/cpsp_day7_vs_spp1_ko_cpsp.cellchat.2.pdf",width=12,height = 6)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht1, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(ht2, newpage = FALSE)
upViewport()
dev.off()

p1<-netAnalysis_contribution(object.list[[1]], signaling = pathways.show)+ggtitle(names(object.list[1]))
p2<-netAnalysis_contribution(object.list[[2]], signaling = pathways.show)+ggtitle(names(object.list[2]))
pdf("figs/cpsp_day7_vs_spp1_ko_cpsp.cellchat.3.pdf",width=12,height = 6)
p1|p2
dev.off()

p1<-netVisual_bubble(object.list[[1]], sources.use = 4, targets.use = c(1:9), remove.isolate = FALSE,signaling =pathways.show,font.size = 20)+ggtitle(names(object.list[1]))
p2<-netVisual_bubble(object.list[[2]], sources.use = 4, targets.use = c(1:9), remove.isolate = FALSE,signaling =pathways.show,font.size = 20)+ggtitle(names(object.list[2]))

pdf("figs/cpsp_day7_vs_spp1_ko_cpsp.cellchat.4.pdf",width=12,height = 6)
p1|p2
dev.off()

plotGeneExpression(cellchat, signaling = pathways.show, split.by = "datasets", colors.ggplot = T, type = "violin")
```
