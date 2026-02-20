import pandas as pd
import anndata as ad

df_expr = pd.read_table('C:\Python\Rawdata_GSE102556\GSE102556_norm_counts_TPM_GRCh38.p13_NCBI.tsv\GSE102556_norm_counts_TPM_GRCh38.p13_NCBI.tsv', sep='\t', index_col=0) 
df_gene_map = pd.read_table('C:\Python\Rawdata_GSE102556\Human.GRCh38.p13.annot.tsv\Human.GRCh38.p13.annot.tsv', sep='\t', index_col=0)
df_sample_info = pd.read_csv('C:\Python\Rawdata_GSE102556\metadata.csv', index_col=0)


df_expr.index = df_expr.index.map(df_gene_map['Symbol'])
df_expr_T = df_expr.T

ofc_samples = df_sample_info[df_sample_info['Tissue'] == 'OFC'].index
df_expr_ofc = df_expr_T.loc[df_expr_T.index.intersection(ofc_samples)]
df_metadata_ofc = df_sample_info.loc[df_expr_ofc.index]


def save_tempo_data(expr_df, meta_df, group_name):
    target_indices = meta_df[meta_df['Phenotype'] == group_name].index
    sub_expr = expr_df.loc[target_indices]
    sub_meta = meta_df.loc[target_indices]
    
    adata = ad.AnnData(X=sub_expr.values, 
                       obs=sub_meta, 
                       var=pd.DataFrame(index=sub_expr.columns))

    if not adata.var_names.is_unique:
        adata.var_names_make_unique()
    
    
    adata.var_names = [name if name != "BMAL1" else "ARNTL" for name in adata.var_names]

    
    adata.write_h5ad(f'C:/Python/Tempo_predata/AnnData/GSE102556/OFC_{group_name}_tempo_input.h5ad')
    print(f"OFC {group_name}: {len(sub_meta)} samples")


save_tempo_data(df_expr_ofc, df_metadata_ofc, 'CTRL')
save_tempo_data(df_expr_ofc, df_metadata_ofc, 'MDD')