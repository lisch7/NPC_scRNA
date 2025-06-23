import omicverse as ov
import scanpy as sc


def clean_adata(adata):
    # Convert adata.raw to an AnnData object
    adata_raw = adata.raw.to_adata()
    
    # Retrieve layers 'counts'
    ov.utils.retrieve_layers(adata_raw, layers='counts')

    # Determine which columns to delete in adata.var
    columns_to_keep = set(['gene_ids', 'feature_types', 'mt', 'ribo', 'hb', 'n_cells', 'percent_cells'])
    all_columns = set(adata_raw.var.columns)
    var_keys_to_delete = list(all_columns - columns_to_keep)

    # Delete specified var columns
    for key in var_keys_to_delete:
        del adata_raw.var[key] 
    
    # Clear uns, obsm, and obsp
    adata_raw.uns.clear()
    adata_raw.obsm.clear()
    adata_raw.obsp.clear()

    # Copy the cleaned data to a new variable and delete the old one
    adata_cleaned = adata_raw.copy()
    del adata_raw

    return adata_cleaned