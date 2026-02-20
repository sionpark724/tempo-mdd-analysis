import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def prep_results(group_name, folder_path):
    path = f"{folder_path}/tempo_results/opt/cycler_gene_prior_and_posterior.tsv"
    df = pd.read_table(path, sep='\t', index_col='gene')
    
    df = df[['phi_loc', 'A_loc', 'mu_loc', 'Q_prob_loc']].copy()
    
    df['phi_hour'] = df['phi_loc'] * (24 / (2 * np.pi))
    df['phi_hour'] = df['phi_hour'] % 24
    
    df['Group'] = group_name
    return df


df_ctrl = prep_results("CTRL", "C:/Python/Tempo_Result/GSE102556/OFC/260203/CTRL_results")
df_mdd = prep_results("MDD", "C:/Python/Tempo_Result/GSE102556/OFC/260203/MDD_results")

comparison_df = pd.concat([df_ctrl, df_mdd])
comparison_df.reset_index(inplace=True)
print(comparison_df.head())


target_genes = ['CLOCK', 'ARNTL', 'PER2', 'RORC', 'NR1D1', 'NR1D2']
df_plot = comparison_df[comparison_df['gene'].isin(target_genes)].copy()

plt.figure(figsize=(12, 8))

sns.stripplot(data=df_plot, x='phi_hour', y='gene', hue='Group', 
              palette={'CTRL': '#1f77b4', 'MDD': '#ff7f0e'}, # 색상 명시
              size=15, alpha=0.7, jitter=False, dodge=True, linewidth=1, edgecolor='white')

for h in range(0, 25, 4):
    plt.axvline(x=h, color='gray', linestyle='--', alpha=0.2)

plt.xlim(0, 24)
plt.xticks(range(0, 25, 2))
plt.title("Core Clock Genes - Phase Shift", fontsize=15)
plt.xlabel("Peak Time (Hour of Day 0-24h)", fontsize=12)
plt.grid(True, axis='y', linestyle=':', alpha=0.5)
plt.legend(title='Group', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig('C:/Python/Tempo_Result/GSE102556/OFC/260203/figures/Fig1_CoreClockGenesShift.png')


target_genes = ['CLOCK', 'ARNTL', 'PER2', 'RORC', 'NR1D1', 'NR1D2']

fig, axes = plt.subplots(len(target_genes)//2, 2, figsize=(12, 4* len(target_genes)//2))
axes_flat = axes.flatten()
x = np.linspace(0, 24, 100)

for i, gene in enumerate(target_genes):
    ax = axes_flat[i]
    for group in ['CTRL', 'MDD']:
        subset = comparison_df[(comparison_df['gene'] == gene) & (comparison_df['Group'] == group)]
        
        if not subset.empty:
            row = subset.iloc[0]
            phi = row['phi_hour']
            A = row['A_loc']
            mu = row['mu_loc']
            
            y = mu + A * np.cos(2 * np.pi * (x - phi) / 24)
            
            color = '#1f77b4' if group == 'CTRL' else '#ff7f0e'
            ax.plot(x, y, label=f"{group} (A={A:.2f})", color=color, linewidth=2)
            ax.fill_between(x, mu, y, color=color, alpha=0.1)
    
    ax.set_title(f"{gene}", fontsize=12)
    ax.set_ylabel("Estimated Expression")
    ax.set_xticks(range(0, 25, 4))
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, linestyle=':', alpha=0.5)

plt.xlabel("Time (Hours)")
plt.tight_layout()
plt.subplots_adjust(hspace=0.6)
plt.savefig('C:/Python/Tempo_Result/GSE102556/OFC/260203/figures/Fig2_CoreClockGenesExpression.png')




comparison_df['amp_diff'] = comparison_df.apply(
    lambda x: comparison_df[(comparison_df['gene'] == x['gene']) & (comparison_df['Group'] == 'CTRL')]['A_loc'].values[0] - 
              comparison_df[(comparison_df['gene'] == x['gene']) & (comparison_df['Group'] == 'MDD')]['A_loc'].values[0] 
    if len(comparison_df[(comparison_df['gene'] == x['gene']) & (comparison_df['Group'] == 'MDD')]) > 0 else 0, axis=1
)


dampening_top = comparison_df[comparison_df['Group'] == 'CTRL'].sort_values(by='amp_diff', ascending=False)
print("\n Dampening genes Top 10:")
print(dampening_top[['gene', 'amp_diff']].head(10))


vanishing_genes = comparison_df[(comparison_df['Group'] == 'MDD') & (comparison_df['A_loc'] < 0.05)]
print("\n Vanishing genes:")
print(vanishing_genes[['gene', 'A_loc']])



plt.figure(figsize=(12, 6))
top_10_genes = dampening_top.head(10)['gene'].tolist()
plot_data = comparison_df[comparison_df['gene'].isin(top_10_genes)]

sns.barplot(data=plot_data, x='gene', y='A_loc', hue='Group', 
            palette={'CTRL': '#1f77b4', 'MDD': '#ff7f0e'}, order=top_10_genes)

plt.title('Rhythm Dampening Genes', fontsize=15)
plt.ylabel('Amplitude', fontsize=12)
plt.xlabel('Gene', fontsize=12)
plt.legend(title='Group', loc='upper right')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('C:/Python/Tempo_Result/GSE102556/OFC/260203/figures/Fig3_DampeningGenes10.png')


plt.figure(figsize=(14, 6))
vanishing_list = vanishing_genes['gene'].tolist()
vanishing_plot_data = comparison_df[comparison_df['gene'].isin(vanishing_list)]

sns.barplot(data=vanishing_plot_data, x='gene', y='A_loc', hue='Group',
            palette={'CTRL': '#1f77b4', 'MDD': '#ff7f0e'})

plt.axhline(y=0.05, color='red', linestyle='--', label='Vanishing Threshold')
plt.title('Rhythm Vanishing Genes', fontsize=15)
plt.ylabel('Amplitude', fontsize=12)
plt.xticks(rotation=45)
plt.legend(title='Group')
plt.tight_layout()
plt.savefig('C:/Python/Tempo_Result/GSE102556/OFC/260203/figures/Fig4_VanishingGenes.png')



comparison_df['phase_diff'] = comparison_df.apply(
    lambda x: comparison_df[(comparison_df['gene'] == x['gene']) & (comparison_df['Group'] == 'MDD')]['phi_hour'].values[0] - 
              comparison_df[(comparison_df['gene'] == x['gene']) & (comparison_df['Group'] == 'CTRL')]['phi_hour'].values[0]
    if len(comparison_df[(comparison_df['gene'] == x['gene']) & (comparison_df['Group'] == 'MDD')]) > 0 else 0, axis=1
)


comparison_df['phase_diff'] = (comparison_df['phase_diff'] + 12) % 24 - 12


significant_shift = comparison_df[(comparison_df['Group'] == 'CTRL') & (comparison_df['phase_diff'].abs() >= 1.5)]
significant_shift = significant_shift.sort_values(by='phase_diff')

print("\n Significant phase shift genes:")
print(significant_shift[['gene', 'phase_diff']])




plt.figure(figsize=(10, 8))
shift_genes = significant_shift['gene'].tolist()
plot_data = comparison_df[comparison_df['gene'].isin(shift_genes)]


for gene in shift_genes:
    gene_data = plot_data[plot_data['gene'] == gene]
    plt.plot(gene_data['phi_hour'], [gene]*2, color='gray', linestyle='--', alpha=0.5, zorder=1)


sns.scatterplot(data=plot_data, x='phi_hour', y='gene', hue='Group', 
                palette={'CTRL': '#1f77b4', 'MDD': '#ff7f0e'}, s=150, zorder=2)

plt.title('Significant Phase Shift Genes', fontsize=15)
plt.xlabel('Peak Time (Hour of Day 0-24h)', fontsize=12)
plt.xticks(range(0, 25, 2))
plt.grid(axis='x', linestyle=':', alpha=0.6)
plt.tight_layout()
plt.savefig('C:/Python/Tempo_Result/GSE102556/OFC/260203/figures/Fig5_PhaseShiftGenes.png')




plt.figure(figsize=(12, 10))

plot_base = comparison_df[comparison_df['Group'] == 'CTRL']

scatter = sns.scatterplot(
    data=plot_base, x='phase_diff', y='amp_diff', 
    hue='Q_prob_loc', size='A_loc', sizes=(50, 500), 
    palette='viridis', alpha=0.7
)

for i, row in plot_base.iterrows():
    if abs(row['phase_diff']) >= 1.5 or row['amp_diff'] >= 0.1:
        plt.text(row['phase_diff'] + 0.1, row['amp_diff'] + 0.005, 
                 row['gene'], fontsize=10, weight='bold')


plt.axvline(x=0, color='black', linestyle='-', alpha=0.3)
plt.axhline(y=0, color='black', linestyle='-', alpha=0.3)
plt.title('Circadian Genes - OFC', fontsize=15)
plt.xlabel('Phase Shift', fontsize=12)
plt.ylabel('Amplitude Dampening', fontsize=12)
plt.grid(True, linestyle=':', alpha=0.4)

plt.savefig('C:/Python/Tempo_Result/GSE102556/OFC/260203/figures/Fig6_Final.png')

comparison_df.to_csv('C:/Python/Tempo_Result/GSE102556/OFC/260203/figures/Comparison_DF')