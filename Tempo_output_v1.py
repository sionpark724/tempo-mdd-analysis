import pandas as pd
import numpy as np


cell_posterior_path = r'C:\Python\Tempo_Result\GSE102556\OFC\260203\MDD_results\tempo_results\opt\cell_posterior.tsv'
cell_dist_df = pd.read_table(cell_posterior_path, sep='\t', index_col='barcode')

print(cell_dist_df.head())


gene_dist_path = r'C:\Python\Tempo_Result\GSE102556\OFC\260203\MDD_results\tempo_results\opt\cycler_gene_prior_and_posterior.tsv'
gene_dist_df = pd.read_table(gene_dist_path, sep='\t', index_col='gene')

print(gene_dist_df.head())


bins = np.linspace(0, 2*np.pi, 24, endpoint=False)


cell_dist_df['theta_max'] = bins[np.argmax(cell_dist_df.iloc[:, :24].values, axis=1)]
print(cell_dist_df[['theta_max']].head())

# core results
core_results = gene_dist_df[['phi_loc', 'A_loc', 'mu_loc', 'Q_prob_loc']]
print(core_results)