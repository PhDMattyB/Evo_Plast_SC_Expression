setwd('~/Parsons_Postdoc/NCC_NEOF_Project/Trailmaker/')

library(Seurat)
library(tidyverse)

sc_data = readRDS(file = 'SSC_Seurat_obj_processed_matrix.rds')

sc_data = readRDS(file = 'Evo_Plast_sc_data.rds')

sc_data@misc
sc_data@meta.data
sc_data$

nebula_data = scToNeb(obj = sc_data, 
                      assay = 'RNA', 
                      id = 'samples', 
                      pred = c('Temp', 
                               'Ecotype', 
                               'Ecotype_temp', 
                               'seurat_clusters', 
                               'cells_id'), 
                      offset = 'nCount_RNA')
