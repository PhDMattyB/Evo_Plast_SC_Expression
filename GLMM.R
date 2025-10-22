setwd('~/Parsons_Postdoc/NCC_NEOF_Project/Trailmaker/')

library(Seurat)
library(tidyverse)
library(nebula)

sc_data = readRDS(file = 'Evo_Plast_sc_data.rds')

sc_data@misc
sc_data@meta.data

umap = DimPlot(sc_data, 
               reduction = 'umap', 
               label = F)

VizDimLoadings(sc_data, dims = 1:2, reduction = 'umap')

sc_data@reductions$umap
sc_data@reductions
sc_data@meta.data %>% 
  as_tibble() %>% 
  write_tsv('SingleCell_Umap_metadata.txt')
sc_data@assays$RNA

sc_data@active.ident

sc_data[['umap']]@feature.loadings
sc_data[['umap']]@key

sc_data@reductions
sc_data[['umap']]@cell.embeddings %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  rowid_to_column() %>% 
  write_csv('SingleCell_Umap_Projections.csv')

sc_data2 = RunPCA(sc_data)
sc_data2 = FindNeighbors(sc_data2)
sc_data2 = FindClusters(sc_data2)
# clust_markers = FindAllMarkers(neighbors,
#                                only.pos = T, 
#                                resolution = 0.5)

ElbowPlot(sc_data2, ndims = 50)
sc_data2 = JackStraw(sc_data2, num.replicate = 100)
sc_data2 = ScoreJackStraw(sc_data2, dims = 1:50)
JackStrawPlot(sc_data2, dims = 1:50)

sc_data2 = RunUMAP(sc_data2, dims = 1:50)

clust_markers = FindAllMarkers(sc_data2,
                               only.pos = T,
                               resolution = 0.5, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)


clust_markers %>%
  # rownames_to_column() %>%
  as_tibble() %>%
#   filter(cluster == 'Cluster 29')%>%
  filter(p_val_adj <= 0.05) %>%
  # arrange(desc(avg_log2FC)) %>%
#   select(gene) %>%
#   slice(1:1000) %>%
  write_tsv('SingleCell_Clusters_All_Markers_UMAP.txt')

FeaturePlot(sc_data2, 
            features = c('pks1'))


DimPlot(sc_data2, 
        reduction = 'umap', 
        label = F)

# Nebula model ------------------------------------------------------------


 nebula_data <- scToNeb(obj = sc_data, 
                      assay = 'RNA', 
                      id = 'samples', 
                      pred = c('Temp', 
                               'Ecotype', 
                               'Ecotype_temp', 
                               'seurat_clusters', 
                               'cells_id'), 
                      offset = 'nCount_RNA')
 
df = model.matrix(~seurat_clusters + Temp + Ecotype + seurat_clusters*Temp + seurat_clusters*Ecotype + Temp*Ecotype + seurat_clusters*Temp*Ecotype,
                  data = nebula_data$pred)

head(df)
dim(df)

# data_grouped = group_cell(count = nebula_data$count, 
#                           id = nebula_data$id, 
#                           pred = df)

re = nebula(nebula_data$count, 
            nebula_data$id, 
            pred = df, 
            ncore = 8)
   