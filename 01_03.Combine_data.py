#!/usr/bin/env python
# coding: utf-8

import configparser
import json
import os

import anndata
import omicverse as ov
import pandas as pd
import scanpy as sc

config = configparser.ConfigParser()
config.read('00-01.config_filtered.ini')

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')
ov.utils.ov_plot_set()

# Defining file and directory paths
merge_data_file = config['Paths']['merge_data_file'] 
merge_data_dir = os.path.dirname(merge_data_file)
os.makedirs(merge_data_dir, exist_ok=True)

filter_data_path = config['Paths']['filter_data_path']


# Ensure that the location where the filtered files are saved exists
if not os.path.exists(filter_data_path):
    raise FileNotFoundError(f"The location to save filtered files {filter_data_path} does not exist.")

# Get the folder name of all samples
try:
    sample_h5ad_files = [os.path.join(filter_data_path, filename) for filename in os.listdir(filter_data_path) if filename.endswith('.h5ad')]
except Exception as e:
    raise Exception(f"Error reading sample folder names: {e}")

# Create a list to store AnnData objects for each sample
adata_list = []
for sample_dir in sample_h5ad_files:
    try:
        adata_list.append(sc.read_h5ad(sample_dir))
    except Exception as e:
        raise Exception(f"Failed to read sample file {sample_dir}: {e}")

concatenated_adata = anndata.concat(adata_list, merge="same")
concatenated_adata.obs_names_make_unique()
concatenated_adata.var_names_make_unique()

def print_size_in_MB(x):
    print('{:.3} MB'.format(x.__sizeof__() / 1e6))

print_size_in_MB(concatenated_adata)

adata = concatenated_adata.copy()
del concatenated_adata

def add_group_info(adata, excel_file, obs_columns, combine_colums):
    try:
        excel_data = pd.read_excel(io=excel_file)
    except Exception as e:
        raise Exception(f"Failed to read Excel file {excel_file}: {e}")

    if combine_colums not in excel_data.columns:
        raise ValueError(f"Column {combine_colums} does not exist in the Excel file.")

    for column in obs_columns:
        if column not in excel_data.columns:
            raise ValueError(f"Column {column} does not exist in the Excel file.")
        group_dict = dict(zip(excel_data[combine_colums], excel_data[column]))
        adata.obs[column] = [group_dict.get(sample, 'NA') for sample in adata.obs['sample']]

excel_store_group = config['add_group']['file_path']
combine_colums = config['add_group']['index_colum']
group_name_str = config['add_group']['group_name']
group_name = json.loads(group_name_str)


# Adding group information
add_group_info(adata, excel_file=excel_store_group, obs_columns=group_name,combine_colums=combine_colums)

# Checking group information
obs_df = adata.obs
selected_columns = obs_df.loc[:, ['sample'] + group_name]
selected_columns = selected_columns.drop_duplicates().reset_index(drop=True)
print(selected_columns)

# Saving selected columns as a CSV file
group_check_path = config['add_group']['group_check_path']
selected_columns.to_csv(f'{group_check_path}/selected_columns.csv', index=False)

# Saving the result file
adata.write_h5ad(merge_data_file)
