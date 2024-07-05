"""
Outputs significant gene expression analysis 
"""
from cross_disorder_copy import cross_disorder_effect_z
from gene_expression import corr_gene_exp_fdr, max_cor_genes
from enigmatoolbox.datasets import fetch_ahba
from enrichr import enrichr

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
print(pca_comp_cor)
umap_comp = components['umap']
umap_comp_cor = umap_comp[2]
print(umap_comp)
r_gene_pca, p_gene_pca = corr_gene_exp_fdr(pca_comp_cor['cortex'][:,0])
r_gene_umap, p_gene_umap = corr_gene_exp_fdr(umap_comp_cor['cortex'][:,0])

# Create dictionary with gene names and correlations
gene_pca_list = [
    {'name':genelabels,'r_gene':r_gene_pca,'p_gene':p_gene_pca}
    for genelabels,r_gene_pca,p_gene_pca in zip(genelabels,r_gene_pca,p_gene_pca)]

gene_umap_list = [
    {'name':genelabels,'r_gene':r_gene_umap,'p_gene':p_gene_umap}
    for genelabels,r_gene_umap,p_gene_umap in zip(genelabels,r_gene_umap,p_gene_umap)]

gene_pca_pos = max_cor_genes(gene_pca_list,0.5,0.05)
gene_pca_neg = max_cor_genes(gene_pca_list,-0.5,0.05)

gene_umap_pos = max_cor_genes(gene_umap_list,0.5,0.05)
gene_umap_neg = max_cor_genes(gene_umap_list,-0.5,0.05)

neg_combined = list(set(gene_umap_neg).intersection(gene_pca_neg))
pos_combined = list(set(gene_umap_pos).intersection(gene_pca_pos))

# enrichr gene analysis
enrichr(neg_combined,"significant negative genes")
enrichr(pos_combined,"signficant postive genes")
