# Code used in the study (Figure 4 & Extended Data Figure 6): "Molecular characterization of colorectal cancer related peritoneal metastatic disease", Kristiaan J. Lenos et al.
# Please direct any questions or comments related to this code to Leandro Moreno: l.ferreiramoreno@amsterdamumc.nl

#Fig 4-H
sample5499 <- Read10X(data.dir = "filtered_feature_bc_matrix_5499/")
sample5502 <- Read10X(data.dir = "filtered_feature_bc_matrix_5502/")
sample5503 <- Read10X(data.dir = "filtered_feature_bc_matrix_5503/")
sample5506 <- Read10X(data.dir = "filtered_feature_bc_matrix_5506/")
sample5504 <- Read10X(data.dir = "filtered_feature_bc_matrix_5504/")
sample5505 <- Read10X(data.dir = "filtered_feature_bc_matrix_5505/")

sample5504_object = CreateSeuratObject(counts = sample5504)
sample5505_object = CreateSeuratObject(counts = sample5505)
sample5499_object = CreateSeuratObject(counts = sample5499)
sample5502_object = CreateSeuratObject(counts = sample5502)
sample5503_object = CreateSeuratObject(counts = sample5503)
sample5506_object = CreateSeuratObject(counts = sample5506)

rm(sample5504, sample5505, sample5499, sample5502, sample5503, sample5506) 

sample5500_object[["percent.mt"]] <- PercentageFeatureSet(sample5500_object, pattern = "^MT-")
sample5501_object[["percent.mt"]] <- PercentageFeatureSet(sample5501_object, pattern = "^MT-")
sample5504_object[["percent.mt"]] <- PercentageFeatureSet(sample5504_object, pattern = "^MT-")
sample5505_object[["percent.mt"]] <- PercentageFeatureSet(sample5505_object, pattern = "^MT-")
sample5499_object[["percent.mt"]] <- PercentageFeatureSet(sample5499_object, pattern = "^MT-")
sample5502_object[["percent.mt"]] <- PercentageFeatureSet(sample5502_object, pattern = "^MT-")
sample5503_object[["percent.mt"]] <- PercentageFeatureSet(sample5503_object, pattern = "^MT-")
sample5506_object[["percent.mt"]] <- PercentageFeatureSet(sample5506_object, pattern = "^MT-")

sample5504_object[["loc"]] <- "primary"
sample5505_object[["loc"]] <- "metastasis"
sample5499_object[["loc"]] <- "metastasis"
sample5502_object[["loc"]] <- "metastasis"
sample5503_object[["loc"]] <- "metastasis"
sample5506_object[["loc"]] <- "metastasis"
sample5504_object[["pat_id"]] <- "B"
sample5505_object[["pat_id"]] <- "B"
sample5499_object[["pat_id"]] <- "C"
sample5502_object[["pat_id"]] <- "D"
sample5503_object[["pat_id"]] <- "E"
sample5506_object[["pat_id"]] <- "F"

sample5504_object_subset <- subset(sample5504_object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
sample5505_object_subset <- subset(sample5505_object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
sample5499_object_subset <- subset(sample5499_object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
sample5502_object_subset <- subset(sample5502_object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
sample5503_object_subset <- subset(sample5503_object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
sample5506_object_subset <- subset(sample5506_object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)

sample5504_object_subset= SCTransform(object = sample5504_object_subset,vars.to.regress = c("nFeature_RNA", "percent.mt"))
sample5505_object_subset= SCTransform(object = sample5505_object_subset,vars.to.regress = c("nFeature_RNA", "percent.mt"))
sample5499_object_subset= SCTransform(object = sample5499_object_subset,vars.to.regress = c("nFeature_RNA", "percent.mt"))
sample5502_object_subset= SCTransform(object = sample5502_object_subset,vars.to.regress = c("nFeature_RNA", "percent.mt"))
sample5503_object_subset= SCTransform(object = sample5503_object_subset,vars.to.regress = c("nFeature_RNA", "percent.mt"))
sample5506_object_subset= SCTransform(object = sample5506_object_subset,vars.to.regress = c("nFeature_RNA", "percent.mt"))


sample5504_object_subset <- RunPCA(sample5504_object_subset, verbose = FALSE)
sample5505_object_subset <- RunPCA(sample5505_object_subset, verbose = FALSE)
sample5504_object_subset <- RunUMAP(sample5504_object_subset, dims = 1:30, verbose = FALSE)
sample5505_object_subset <- RunUMAP(sample5505_object_subset, dims = 1:30, verbose = FALSE)

sample5504_object_subset <- FindNeighbors(sample5504_object_subset, dims = 1:30, verbose = FALSE)
sample5505_object_subset <- FindNeighbors(sample5505_object_subset, dims = 1:30, verbose = FALSE)

sample5504_object_subset <- FindClusters(sample5504_object_subset, verbose = FALSE, resolution = 0.6)
sample5505_object_subset <- FindClusters(sample5505_object_subset, verbose = FALSE, resolution = 0.6)

#Fig 4-G
merged_objects_all=merge(x = sample5504_object,y = c(sample5505_object, sample5499_object, sample5502_object, sample5503_object, sample5506_object))
merged_objects_all= SCTransform(object = merged_objects_all,vars.to.regress = c("nFeature_RNA"))

merged_objects_all <- RunPCA(merged_objects_all, verbose = FALSE)
merged_objects_all <- RunUMAP(merged_objects_all, dims = 1:10, verbose = FALSE)
merged_objects_all <- FindNeighbors(merged_objects_all, dims = 1:10, verbose = FALSE)
merged_objects_all <- FindClusters(merged_objects_all, verbose = FALSE, resolution = 0.2)
umap_celcius_allpatient <- DimPlot(merged_objects_all, label = FALSE,  group.by="pat_id", cols = c('B-Primary' = "#ff75a0", "B-Metastasis" = "#93329e",'C' = "#f0a500", 'D' = "#1687a7", 'E' = "#666666", "F" = "#95e1d3"))

#Fig 4-H
merged_objects=merge(x = sample5504_object_subset,y = sample5505_object_subset)
merged_objects= SCTransform(object = merged_objects,vars.to.regress = c("nFeature_RNA"))

merged_objects <- RunPCA(merged_objects, verbose = FALSE)
merged_objects <- RunUMAP(merged_objects, dims = 1:10, verbose = FALSE, metric = "euclidean")
merged_objects <- FindNeighbors(merged_objects, dims = 1:10, verbose = FALSE)
merged_objects <- FindClusters(merged_objects, verbose = FALSE, resolution = 0.05)
umap_celcius <- DimPlot(merged_objects, label = TRUE, cols = c('0' = "#e5707e", '1' = "#e6b566", '2' = "#e8e9a1", '3' = "#a3ddcb", "4" = "#16c786",  "5" = "#1687a7"))
umap_loc <- DimPlot(merged_objects, label = FALSE,  group.by="loc", cols = c("primary" = "#ff75a0", "metastasis" = "#00917c"))  

ggsave("umap_loc_patB.pdf", plot = umap_loc, dpi=300, width=10, height=8)


##Extended data Fig 6-A,B
#Extended Fig 6-A; adaptated from StackedVlnPlot

library(patchwork)

modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)),
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


markers <- StackedVlnPlot(obj = merged_objects_all, idents = c("Stromal cells","Myeloid cells","T cells","B cells","Endothelial cells"),
                               features = c('COL3A1','COL1A2','COL1A1', 'HIF1A','GLUL','PLXDC2', 'BANK1', 'ST6GAL1', 'FCHSD2', 'FYN', "PTPRC", "PRKCH", 'CALCRL', 'GRB10', 'HSPG2' ), 
                               cols =  c('Stromal cells' = "#b82969", 'Myeloid cells' = "#d4b920", 'T cells' = "#226589", 'B cells' = "#636d7a", "Endothelial cells" = "#16c786"))

ggsave("C:/amc.intra/users/L/lferreiramoreno/home/R_plots/markers.pdf", plot = markers_peri, dpi=300, width=10, height=20)


#Extended Fig 6-B
library(dplyr)
merged_objects.markers <- FindAllMarkers(merged_objects, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
merged_objects.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC )
top10 <- merged_objects.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC )

library(viridis)
heatmap_patienta <- DoHeatmap(merged_objects, features = top10$gene, size=3.5, group.colors = c('0' = "#e5707e", '1' = "#e6b566", '2' = "#e8e9a1", '3' = "#a3ddcb", "4" = "#16c786", "5" = "#1687a7")) +  theme(text = element_text(size = 10)) + scale_fill_viridis()
ggsave("C:/amc.intra/users/L/lferreiramoreno/home/R_plots/heatmap_top10_patB.pdf", plot = heatmap_patienta,
       dpi=300, width=12, height=10)




