setwd('~/Parsons_Postdoc/NCC_NEOF_Project/Trailmaker/')

library(Seurat)
library(tidyverse)
library(nebula)

sc_data = readRDS(file = 'Evo_Plast_sc_data.rds')

sc_data@misc
sc_data@meta.data


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

