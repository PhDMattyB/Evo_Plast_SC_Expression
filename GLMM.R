setwd('~/Parsons_Postdoc/NCC_NEOF_Project/Trailmaker/')

library(Seurat)
library(tidyverse)
library(nebula)

sc_data = readRDS(file = 'Evo_Plast_sc_data.rds')

sc_data@misc
sc_data@meta.data

umap = DimPlot(sc_data, 
               reduction = 'umap', 
               label = T)

clust_markers = FindAllMarkers(sc_data, 
                               only.pos = T)

clust_markers %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  filter(cluster == 'Cluster 29')%>%
  filter(p_val_adj <= 0.05) %>%
  arrange(desc(avg_log2FC)) %>% 
  select(gene) %>%
  slice(1:1000) %>%
  write_tsv('Cluster29_T1000_Marker_genes.txt')

 nebula_data <- scToNeb(obj = sc_data, 
                      assay = 'RNA', 
                      id = 'samples', 
                      pred = c('Temp', 
                               'Ecotype', 
                               'Ecotype_temp', 
                               'seurat_clusters', 
                               'cells_id'), 
                      offset = 'nCount_RNA')
 
# df = model.matrix(~seurat_clusters+cells_id, 
#                   data = nebula_data$pred)

