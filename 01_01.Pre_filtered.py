#!/usr/bin/env python
# coding: utf-8

import configparser
import os
import sys
from distutils.util import strtobool

import omicverse as ov
import scanpy as sc


config = configparser.ConfigParser()
config.read('00-01.config_filtered.ini')

# 1.2 scanpy setting #
sc.settings.verbosity = 3   
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, facecolor='white')

#### 2.1 set the reading path ####

data_dir = config['Paths']['raw_data_path'] 
sample_name =  sys.argv[1] # import the sample name

sample_data_path = os.path.join(data_dir,sample_name,'outs', 'filtered_feature_bc_matrix')

adata = sc.read_10x_mtx(sample_data_path,  ## the directory with the `.mtx` file
                var_names='gene_symbols', ## use gene symbols for the variable names (variables-axis index)
                cache=True ## if True, write a cache file for faster subsequent reading
                       )

adata.obs_names=[sample_name+ '_' +x for x in adata.obs_names]

#### filter low quality cells ========================================================================
#### Add the sample name as a new obs
adata.obs['sample'] = sample_name


# sc.pl.highest_expr_genes(adata, n_top=20, )

#### mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")

#### ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

#### hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )

#### cell_complexity
# adata.obs['cell_complexity'] = adata.obs['n_genes_by_counts'].map(math.log10) / adata.obs['total_counts'].map(math.log10)
# adata

minGene = int(config['filter']['minGene'])
maxGene = int(config['filter']['maxGene'])
minCount = int(config['filter']['minCount'])
maxCount = int(config['filter']['maxCount'])
pct_MT = float(config['filter']['pct_MT'])
pct_HB = float(config['filter']['pct_HB']) 

adata = adata[
    (adata.obs.n_genes_by_counts >= minGene) &
    (adata.obs.n_genes_by_counts <= maxGene) &
    (adata.obs.total_counts >= minCount) &
    (adata.obs.total_counts <= maxCount) &
    (adata.obs.pct_counts_mt < pct_MT) &
    (adata.obs.pct_counts_hb < pct_HB),
    :]

#### Enviroment Correct ==========================================================================================
Env_correct=config.getboolean('env_correct', 'env_correct_used')
Env_method=config['env_correct']['env_correct_method']

import logging

import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import scipy.sparse
from rpy2.robjects import FloatVector, IntVector, r
from rpy2.robjects.packages import importr

rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()

if Env_correct:
    if Env_method == "Soupx":
        adata_pp = adata.copy()
        sc.pp.normalize_per_cell(adata_pp)
        sc.pp.log1p(adata_pp)
        sc.pp.pca(adata_pp)
        sc.pp.neighbors(adata_pp)
        sc.tl.leiden(adata_pp, key_added="soupx_groups")
        soupx_groups = adata_pp.obs["soupx_groups"]
        del adata_pp
        cells = adata.obs_names
        genes = adata.var_names
        data = adata.X.T
        sample_raw_martix_path = os.path.join(data_dir,sample_name,'outs', 'raw_feature_bc_matrix')
        adata_raw = sc.read_10x_mtx(sample_raw_martix_path,var_names='gene_symbols',cache=True)
        adata_raw.var_names_make_unique()
        adata_raw.obs_names=[sample_name+ '_' +x for x in adata_raw.obs_names]
        data_tod = adata_raw.X.T

        ro.globalenv['data'] = data
        ro.globalenv['data_tod'] = data_tod
        ro.globalenv['genes'] = genes
        ro.globalenv['cells'] = cells
        ro.globalenv['soupx_groups'] = soupx_groups
 
        ro.r('''
        library(SoupX)

        # specify row and column names of data
        rownames(data) <- genes
        colnames(data) <- cells

        # ensure correct sparse format for table of counts and table of droplets
        data <- as(data, "sparseMatrix")
        data_tod <- as(data_tod, "sparseMatrix")

        # Generate SoupChannel Object for SoupX 
        sc <- SoupChannel(data_tod, data, calcSoupProfile=FALSE)

        # Add extra meta data to the SoupChannel object
        soupProf <- data.frame(row.names=rownames(data), est=rowSums(data)/sum(data), counts=rowSums(data))
        sc <- setSoupProfile(sc, soupProf)
        sc <- setClusters(sc, as.character(soupx_groups))

        # Estimate contamination fraction
        sc <- tryCatch(
        SoupX::autoEstCont(sc),
        error = function(e) {
            print("autoEstCont Error !")
            setContaminationFraction(sc, 0.2)
        }
        )

        # Infer corrected table of counts and rount to integer
        out <- adjustCounts(sc, roundToInt=TRUE)
        ''')

        out = ro.r['out']

        adata.layers["counts"] = adata.X
        adata.layers["soupX_counts"] = out.T
        adata.X = adata.layers["soupX_counts"]

        print(f"Total number of genes: {adata.n_vars}")

        # Min 20 cells - filters out 0 count genes
        sc.pp.filter_genes(adata, min_cells=20)
        print(f"Number of genes after cell filter: {adata.n_vars}")
    else:

        data = adata.X.T

        sample_raw_martix_path = os.path.join(data_dir,sample_name,'outs', 'raw_feature_bc_matrix')
        adata_raw = sc.read_10x_mtx(sample_raw_martix_path,var_names='gene_symbols',cache=True)
        adata_raw.var_names_make_unique()
        adata_raw.obs_names=[sample_name+ '_' +x for x in adata_raw.obs_names]
        data_raw = adata_raw.X.T

        ro.globalenv['data'] = data
        ro.globalenv['data_raw'] = data_raw

        ro.r('''
        library(decontX)
        sce <- SingleCellExperiment(list(counts = data))
        sce.raw <- SingleCellExperiment(list(counts = data_raw))
        sce <- decontX(sce, background = sce.raw)

        out_count = round(decontXcounts(sce))
        decontX_contamination <- sce$decontX_contamination
        ''')

        out_count = ro.r['out_count']
        decontX_contamination = ro.r['decontX_contamination']

        adata.obs["decontX_contamination"] = decontX_contamination

        adata.layers["counts"] = adata.X
        adata.layers["decontX"] = out_count.T
else:
    print("")
    

#### Doublet Detection ============================================================================================
#### scrublet
sc.external.pp.scrublet(adata, random_state= 100)

#### scDblFinder
ro.r('''
suppressMessages(library(Seurat))
suppressMessages(library(scater))
suppressMessages(library(scDblFinder))
suppressMessages(library(BiocParallel))
''')

data_mat = adata.X.T

ro.globalenv['data_mat']=data_mat

ro.r('''
set.seed(123)
sce = scDblFinder(
    SingleCellExperiment(
        list(counts=data_mat),
    ) 
)
doublet_score = sce$scDblFinder.score
doublet_class = sce$scDblFinder.class
''')

doublet_score=ro.r('doublet_score')
doublet_class=ro.r('doublet_class')


adata.obs["scDblFinder_score"] = doublet_score
adata.obs["scDblFinder_class"] = doublet_class
print(adata.obs.scDblFinder_class.value_counts())

#### 5. save the filtered annadata ###########################################################

# set the save path
save_path = config['Paths']['filter_data_path'] 

# Check if the path exists, if not, create it
os.makedirs(save_path, exist_ok=True)

# save the filtered file
results_file = os.path.join(save_path,f"seu_{sample_name}_quality_control.h5ad")
adata.write(results_file)

