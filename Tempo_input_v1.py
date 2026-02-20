import anndata
import pandas as pd
import torch
import numpy as np
import random

# --- Seed ---
seed_value = 2026
random.seed(seed_value)
np.random.seed(seed_value)
torch.manual_seed(seed_value)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(seed_value)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

adata = anndata.read_h5ad('C:\Python\Tempo_predata\AnnData\GSE102556\OFC_MDD_tempo_input.h5ad')


print(adata)
print("--- Sample ---")
print(adata.obs.head())
print("--- Gene ---")
print(adata.var.head())
print("--- Expression Level ---")
print(adata.to_df().iloc[:5, :5])


folder_out = 'C:/Python/Tempo_Result/GSE102556/OFC/260203/MDD_results'


gene_acrophase_prior_path = 'C:\Python\Tempo_predata\core_clock_and_ubiq_acrophase_prior.csv'


core_clock_gene_path = 'C:\Python\Tempo_predata\core_clock_genes.txt'


reference_gene = 'ARNTL' 

min_gene_prop = 0


import tempo
from tempo import unsupervised_alg

tempo.unsupervised_alg.run(
    adata = adata,
    folder_out = folder_out,
    gene_acrophase_prior_path = gene_acrophase_prior_path,
    core_clock_gene_path = core_clock_gene_path,
    reference_gene = reference_gene,
    min_gene_prop = min_gene_prop
)