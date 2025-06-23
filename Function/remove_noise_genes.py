import numpy as np
import pandas as pd


def remove_noise_genes(adata, mito_genes=True, heatshock_genes=True, ribo_genes=True, disso_genes=True,hb_genes=True,TCR_genes=True,BCR_genes=True):
   if mito_genes:
       # Get Mito genes
       mito_data = pd.read_excel("./data_files/gene_list/Mitochondria_genes.xlsx")
       mito_genes_to_exclude = mito_data['gene'].tolist()
   else:
       mito_genes_to_exclude = []
   
   if heatshock_genes:
       # Get Heatshock genes
       heatshock_data = pd.read_excel("./data_files/gene_list/Heat_shock_protein_genes.xlsx")
       heatshock_gene_to_exclude = heatshock_data['gene'].tolist()
   else:
       heatshock_gene_to_exclude = []  
       
   if ribo_genes:
       # Get Ribosome genes
       ribo_data = pd.read_excel("./data_files/gene_list/Ribosome_genes.xlsx")
       ribo_genes_to_exclude = ribo_data['gene'].tolist()
   else:
       ribo_genes_to_exclude = [] 
       
   if disso_genes:
       # Get Dissociation genes
       disso_data = pd.read_excel("./data_files/gene_list/Dissociation_genes.xlsx")
       disso_genes_to_exclude = disso_data['gene'].tolist()
   else:
       disso_genes_to_exclude = []    

   if hb_genes:
       # Get Dissociation genes
       hb_data = pd.read_excel("./data_files/gene_list/Hemoglobin_genes.xlsx")
       hb_genes_to_exclude = hb_data['gene'].tolist()
   else:
       hb_genes_to_exclude = []  

   if TCR_genes:
       # Get Dissociation genes
       TCR_data = pd.read_excel("./data_files/gene_list/TCR_genes.xlsx")
       TCR_genes_to_exclude = TCR_data['gene'].tolist()
   else:
       TCR_genes_to_exclude = [] 
    
   if BCR_genes:
       # Get Dissociation genes
       BCR_data = pd.read_excel("./data_files/gene_list/BCR_genes.xlsx")
       BCR_genes_to_exclude = BCR_data['gene'].tolist()
   else:
       BCR_genes_to_exclude = [] 
                  
   # Combine all genes to exclude
   genes_to_exclude = mito_genes_to_exclude + heatshock_gene_to_exclude + ribo_genes_to_exclude + disso_genes_to_exclude + hb_genes_to_exclude + TCR_genes_to_exclude +  BCR_genes_to_exclude  
   
   genes_to_remove_mask = np.isin(adata.var_names, genes_to_exclude)
   
   adata.var["highly_variable_features"][genes_to_remove_mask] = False