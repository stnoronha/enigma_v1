"""
Outputs significant gene expression analysis 
"""
from cross_disorder_copy import cross_disorder_effect_z
from gene_expression import corr_gene_exp_fdr, max_cor_genes, store_list
from enigmatoolbox.datasets import fetch_ahba,risk_genes

# Get cross disorder effect values from dimensionality reduction methods

dim_reduction = ['pca','umap']
components = {x:cross_disorder_effect_z(method=x) for x in dim_reduction}

genes = fetch_ahba()
reglabels = genes['label']

# remove the label column
genes = genes.drop('label', axis=1)
genelabels = list(genes.columns)

# run gene analysis
pca_comp = components['pca']
pca_comp_cor = pca_comp[0]
umap_comp = components['umap']
umap_comp_cor = umap_comp[2]
r_gene_pca, p_gene_pca = corr_gene_exp_fdr(pca_comp_cor['cortex'][:,0])
r_gene_umap, p_gene_umap = corr_gene_exp_fdr(umap_comp_cor['cortex'][:,0])

# Create dictionary with gene names and correlations
gene_pca_list = [
    {'name':genelabels,'r_gene':r_gene_pca,'p_gene':p_gene_pca}
    for genelabels,r_gene_pca,p_gene_pca in zip(genelabels,r_gene_pca,p_gene_pca)]

gene_umap_list = [
    {'name':genelabels,'r_gene':r_gene_umap,'p_gene':p_gene_umap}
    for genelabels,r_gene_umap,p_gene_umap in zip(genelabels,r_gene_umap,p_gene_umap)]

gene_pca = max_cor_genes(gene_pca_list,0.05)

gene_umap = max_cor_genes(gene_umap_list,0.05)

combined = list(set(gene_pca).intersection(gene_umap))

store_list(gene_pca,'PCA_fdr_significant')

store_list(gene_umap,'UMAP_fdr_significant')

store_list(combined,"combined_fdr_significant")

#compare fdr significant genes to GWAS risk genes for disorders

GWASgenes = {}
for d in ['adhd', 'asd', 'bipolar', 'depression', 'epilepsy', 'hippocampal_volume', 'schizophrenia', 'tourette']:
    GWASgenes[d]=risk_genes(d)

sig_atrisk = {x: set(combined).intersection(GWASgenes[x]) for x in GWASgenes.keys()}

with open('atrisk.txt', 'w') as f:
    for k, v in sig_atrisk.items():
        print(k, v,file=f)
