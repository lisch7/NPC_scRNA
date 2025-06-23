#!/usr/bin/env python
# coding: utf-8

#### 1. preprocess ####

import ast
import configparser
import datetime
import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import omicverse as ov
import pandas as pd
import scanpy as sc

from Function.clean_adata import *
from Function.formatted_title import *
from Function.remove_noise_genes import *

ov.ov_plot_set()

print(f'omicverse version: {ov.__version__}')
print(f'scanpy version: {sc.__version__}')

# Reading the configuration file
config = configparser.ConfigParser()
config.read('00-02.config_preprocess.ini')

print_title_with_time(title="The preprocessing process begins!")

#### Reading settings ############################################################
## Reading and saving files related
output_path = config['Paths']['output_path'] ## Path to save h5ad files
input_file = config['Paths']['input_file'] ## Input file to read
fig_path = config['Paths']['fig_path'] ## Path to save images
save_path = config['Paths']['save_path'] ## Path to save other files
sc.settings.figdir = fig_path

## Normalization settings
hvg=int(config['preprocess']['hvg'])  ## Highly variable genes setting, default 2000
mode_selected=config['preprocess']['mode_selected']  ## Mode for normalization and finding highly variable genes
target_sum_set=int(config['preprocess']['target_sum']) 
no_cc_set = config.getboolean('preprocess', 'no_cc') ## Whether to exclude cell cycle-related genes

batch_key_from_config = config.get('preprocess', 'batch_key')

# Check if batch_key is the string 'None', if so, convert it to Python's None object
if batch_key_from_config.lower() == 'none':
    batch_key = None
else:
    batch_key = batch_key_from_config

## scale settings
regression_status = config.getboolean('scale', 'regression_status') ## Whether to perform regression
regression_condition_str = config['scale']['regression_condition']
regression_condition = ast.literal_eval(regression_condition_str) ## Extract the specific information for regression
n_jobs = config['scale']['n_jobs'] ## Number of cores for scaling

## Dimensionality reduction settings
reduction_method = config['reduction']["method"]
reduction_resolution=float(config['reduction']['reduction_resolution'])
pc_dim = int(config['reduction']['pc_num'])

## Generated file names
mode_selected_name = mode_selected.replace('|', '_')  # Replace | with _
file_name_preprocess = f'seu_preprocess_{mode_selected_name}_hvg{hvg}.h5ad'
file_name_final = f'seu_preprocess_{mode_selected_name}_hvg{hvg}_pc{pc_dim}_{reduction_method}_res{reduction_resolution}.h5ad'

## Generated file paths
file_path_preprocess=os.path.join(output_path, file_name_preprocess) ## File name after normalization and scaling
file_path_final = os.path.join(output_path, file_name_final) ## File name after script completion


s_genes_file = pd.read_excel("./data_files/s_genes.xlsx")
s_genes = s_genes_file['gene'].tolist()

g2m_genes_file = pd.read_excel("./data_files/g2m_genes.xlsx")
g2m_genes = g2m_genes_file['gene'].tolist()

#### Normalization and scaling ########################################################
if os.path.exists(file_path_preprocess):
    print_title_without_time(title="Normalization and scale had finished before!")
else:
    print_title_with_time(title="Normalization and scale begin!")
    
    adata = sc.read_h5ad(input_file)
    ov.utils.store_layers(adata, layers='counts')

    # Normalization and finding highly variable genes
    adata = ov.pp.preprocess(adata, mode=mode_selected, target_sum=target_sum_set, n_HVGs=hvg,no_cc=no_cc_set,batch_key=batch_key)
    adata.raw = adata
    # Calculating cell cycle scores
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes, copy=False,use_raw=True)
    adata = adata[:, adata.var.highly_variable_features]

    # scale   
    if regression_status:
        # Perform regression and scaling
        adata_mock = sc.pp.regress_out(adata, regression_condition, n_jobs=n_jobs, copy=True)
        adata.layers['regressed'] = adata_mock.X
        del adata_mock

        # Scale
        adata_mock = adata.copy()
        adata_mock.X = adata_mock.layers['regressed']
        ov.pp.scale(adata_mock)
        adata.layers['regressed_and_scaled'] = adata_mock.layers['scaled']
        adata.layers['scaled'] = adata_mock.layers['scaled']
        del adata_mock
    else:
        # Only perform scaling
        ov.pp.scale(adata)

    adata.write_h5ad(file_path_preprocess)

print_title_with_time(title="Normalization and scale completed!")        

#### 3. Dimensionality Reduction ####
print_title_with_time(title="Reduction begin")  

## Generated storage names
reduction_store_name = f"{reduction_method}_res{reduction_resolution}"

if os.path.exists(file_path_final):
    print_title_without_time(title="Reduction has finished before!") 
else:
    if 'adata' in globals():
        print("adata already exists in the environment. Skipping read_h5ad!")
    else:
        adata = sc.read_h5ad(file_path_preprocess)
    
    ov.pp.pca(adata, layer='scaled', n_pcs=50)    
    ov.utils.plot_pca_variance_ratio(adata)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=pc_dim, use_rep='scaled|original|X_pca',key_added="neighbors_original")
    sc.tl.umap(adata,neighbors_key="neighbors_original")
    adata.obsm['X_umap_original']=adata.obsm['X_umap']
    ov.utils.cluster(adata,method=reduction_method,key_added=reduction_store_name, neighbors_key="neighbors_original",resolution=reduction_resolution)
    adata.write_h5ad(file_path_final)
    
print_title_with_time(title="Reduction completed and final file saved successfully!") 


#### 4. Plotting dimensionality reduction and gene plots ####
print_title_with_time(title="Drawing the required pictures") 

with plt.rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata,color=[reduction_store_name],neighbors_key="neighbors_harmony",frameon=False, legend_loc='on data',
               legend_fontsize=12, legend_fontoutline=2)

# UMAP and Featureplot
umap_name = 'X_umap_original'
ov.utils.embedding(adata,basis= umap_name,color={reduction_store_name},frameon='small',show=False,save=f"_{reduction_store_name}.pdf") # UMAP图
ov.utils.embedding(adata,basis= umap_name,color=["CD3D","CD3E","CD3G","CD4","FOXP3"],frameon='small',show=False,save=f"_{reduction_store_name}_CD4T.pdf") # CD4T-cell
ov.utils.embedding(adata,basis= umap_name,color=["CD3D","CD3E","CD3G","CD8A","KLRD1","KLRF1"],frameon='small',show=False,save=f"_{reduction_store_name}_CD8T.pdf") # CD8T-cell
ov.utils.embedding(adata,basis= umap_name,color=["CD3D","CD3E","CD3G","CD8A","KLRD1","KLRF1","NKG7"],frameon='small',show=False,save=f"_{reduction_store_name}_T_cell.pdf") # NK-cell
ov.utils.embedding(adata,basis= umap_name,color=["CD79A","CD79B","CD19","MS4A1","MZB1"],frameon='small',show=False,save=f"_{reduction_store_name}_B_cell.pdf") # B-cell
ov.utils.embedding(adata,basis= umap_name,color=["LYZ","CD68","MS4A6A","CD1E","IL3RA","LAMP3"],frameon='small',show=False,save=f"_{reduction_store_name}_Myeloid.pdf") # Myeloid
ov.utils.embedding(adata,basis= umap_name,color=['CSF1R','CSF3R','S100A8','S100A9'],frameon='small',show=False,save=f"_{reduction_store_name}_Neutrophils.pdf") # Neu
ov.utils.embedding(adata,basis= umap_name,color=["CLDN5","VWF"],frameon='small',show=False,save=f"_{reduction_store_name}_Endothelial_cell.pdf") # Endo
ov.utils.embedding(adata,basis= umap_name,color=["COL1A1","COL3A1","DCN"],frameon='small',show=False,save=f"_{reduction_store_name}_Fibroblasts.pdf") # Fibro
ov.utils.embedding(adata,basis= umap_name,color=["EPCAM","KRT19","KRT18","KRT5","KRT15"],frameon='small',show=False,save=f"_{reduction_store_name}_Epithelial_cell.pdf") # Epi
ov.utils.embedding(adata,basis= umap_name,color=["CLEC4C","LILRA4","LILRB4"],frameon='small',show=False,save=f"_{reduction_store_name}_pDC.pdf") # pDC
ov.utils.embedding(adata,basis= umap_name,color=["CPA3","TPSB2","TPSAB1"],frameon='small',show=False,save=f"_{reduction_store_name}_Mast.pdf") # Mast
ov.utils.embedding(adata,basis= umap_name,color=["MKI67","TOP2A","STMN1"],frameon='small',show=False,save=f"_{reduction_store_name}_Proliferating.pdf") # Proliferating
ov.utils.embedding(adata,basis= umap_name,color=["XCR1", "CLEC9A", "BATF3", "IRF8"],frameon='small',show=False,save=f"_{reduction_store_name}_cDC1.pdf") # cDC1
ov.utils.embedding(adata,basis= umap_name,color=["CLEC10A", "FCER1A", "CD1C", "CD1E"],frameon='small',show=False,save=f"_{reduction_store_name}_cDC2.pdf") # cDC2
ov.utils.embedding(adata,basis= umap_name,color=["LAMP3","CCR7"],frameon='small',show=False,save=f"_{reduction_store_name}_cDC3.pdf") # cDC3
ov.utils.embedding(adata,basis= umap_name,color=["CD207", "CD1A", "MAP2K1"],frameon='small',show=False,save=f"_{reduction_store_name}_Langerhans.pdf") # Langerhans
ov.utils.embedding(adata,basis= umap_name,color=["AXL", "SIGLEC6", "PPP1R14A", "KLF4", "ATF5", "IL4I1", "CD5", "CX3CR1", "THBD"],frameon='small',show=False,save=f"_{reduction_store_name}_tDC.pdf") # tDC
ov.utils.embedding(adata,basis='X_umap_original',color=["IGLV3-19","IGHG1","IGKV3-20","HBB"],frameon='small',show=False,save=f"_{reduction_store_name}_polution.pdf") # polution


#### Marker plot
marker_method = config['marker_detection']['marker_detection_method']
marker_method_safe_name = marker_method.replace('-', '_')
marker_key_add = f"{reduction_method}_{marker_method_safe_name}"

sc.tl.dendrogram(adata,reduction_store_name,use_rep='scaled|original|X_pca')
sc.tl.rank_genes_groups(adata, reduction_store_name, use_rep='scaled|original|X_pca',method = marker_method,use_raw=False,key_added=marker_key_add)
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key=marker_key_add, show=False, save=f"_rank_marker_{marker_key_add}.pdf")

#### dotplot and violin plot

marker_genes={
    "T cell":["PTPRC","CD3D","CD4","CD8A"],
    'Naive':['SELL',"CCR7","TCF7","LEF1"],
    'Cytokines_effector':['GZMK','GZMA','GZMB','GZMH','IFNG','NKG7','PRF1','GNLY','IL7R','CD40LG'],
    'Tem':['ANXA1','ANKRD28','CD69'],
    'Resident_marker':['ZNF683','ITGA1','ITGAE'],
    'Co_stimulatory_marker' :['CD28','TNFRSF14',"ICOS","TNFRSF9"],
    'Chemokine':["CXCR5","CXCL13","CXCR3","CXCR6","IL10","CX3CR1"],
    'Immune_checkpoint':["CD69","SLAMF6",'CTLA4','HAVCR2',"BTLA","LAG3","TIGIT","PDCD1","TOX","EOMES","TBX21"],
    'G protein':["GPR183","S1PR1","S1PR5"],
    'Treg':['FOXP3','IL2RA','IKZF2','NCR1'],
    'Proliferative signal':['MKI67','TOP2A','STMN1'],
    'NK':['TYROBP','KLRD1','KLRF1','KLRB1','GNLY','NKG7',"NCAM1"],
    'Tγδ':['TRDV2','TRGV9','TRGC2'],
    'ILC1':['CD200R1','CD5','TMEM176A','TMEM176B'],
    'ILC2':['IL1RL1','HPGDS','IL17RB'],
    'ILC3':['IL1R1','RORC','NRP1','LIF','RUNX2','PCDH9'],
    'Other major celltype':["CD68", "CD163",'MS4A1',"CD79B","MZB1", "JCHAIN" ],
    "IFN":["ISG15","IFIT3","IFIT2"],
    "Other":["FCGR3A","XCL1", "XCL2","SLC4A10","CD160"]
}

adata_raw=adata.raw.to_adata().copy()

marker_genes_in_data = dict()
for ct, markers in marker_genes.items():
    markers_found = list()
    for marker in markers:
        if marker in adata_raw.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found
    
sc.tl.dendrogram(adata_raw,'louvain_res1.0',use_rep="X_pca_harmony")

sc.pl.dotplot(
    adata_raw,
    groupby = reduction_store_name,
    var_names = marker_genes_in_data,
    dendrogram = True,
    standard_scale = "var",  # standard scale: normalize each gene to range from 0 to 1
    save = "majoy_celltype_marker_dotplot.pdf"
)

sc.pl.stacked_violin(adata_raw,
                     groupby = reduction_store_name,
                     var_names = marker_genes_in_data,
                     dendrogram = True,
                     standard_scale = "var",  # standard scale: normalize each gene to range from 0 to 1
                     save = "majoy_celltype_marker_stacked_violin.pdf"
                    )


#### 4. Subpopulation Test ####
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import scipy.sparse
from rpy2.robjects import FloatVector, IntVector, r
from rpy2.robjects.packages import importr

batch_key_reduction = config['reduction']['batch_key_reduction']
cores = int(config['reduction']['parallel_cores'])

