# gene expression analysis
from enigmatoolbox.datasets import fetch_ahba
import matplotlib.pyplot as plt
from enigmatoolbox.utils import parcel_to_surface
from cross_disorder_copy import cross_disorder_effect_z
import numpy as np
from enigmatoolbox.datasets.base import load_summary_stats
from enigmatoolbox.utils.parcellation import parcel_to_surface, surface_to_parcel
from scipy.stats import spearmanr
import csv

# PCA cross disorder effect
components_z, variance_z, names_z = cross_disorder_effect_z()

# UMAP cross disorder effect
umap_cor, umap_sub, umap_comp, names = cross_disorder_effect_z(method="umap")

genes = fetch_ahba()
reglabels = genes['label']

# remove the label column
genes = genes.drop('label', axis=1)
genelabels = list(genes.columns)


def corr_gene_exp(solution):
    if len(solution) == 200:
        # put it on a surface
        d_mri_surf = parcel_to_surface(solution, f'schaefer_200_conte69')
        # map it back on DK parcels used by ENIGMA
        d_mri = surface_to_parcel(d_mri_surf, f'aparc_conte69')
        # remove midbrain
        d_mri = np.delete(d_mri, np.array([0, 4, 39]))
    elif len(solution) == 68:
        d_mri = solution
    else:
        raise ValueError('Solution must be 200 or 68 long')
    r_gene = [] #R is correlation coefficient
    p_gene = [] #p is probability that correlation exists if null is true
    for i, gene in enumerate(genelabels):
        geneexp = genes.iloc[0:68,i]
        # get the location of non-nan values
        idx = ~np.isnan(geneexp)        
        c_gene = spearmanr(d_mri[idx], geneexp[idx])
        r_gene = np.append(r_gene, c_gene[0])
        p_gene = np.append(p_gene, c_gene[1])
    return r_gene, p_gene

r_gene_pca, p_gene_pca = corr_gene_exp(components_z['cortex'][:,0])
r_gene_umap, p_gene_umap = corr_gene_exp(umap_comp['cortex'][:,0])

# Create dictionary with gene names and correlations
gene_pca_list = [
    {'name':genelabels,'r_gene':r_gene_pca,'p_gene':p_gene_pca}
    for genelabels,r_gene_pca,p_gene_pca in zip(genelabels,r_gene_pca,p_gene_pca)]

gene_umap_list = [
    {'name':genelabels,'r_gene':r_gene_umap,'p_gene':p_gene_umap}
    for genelabels,r_gene_umap,p_gene_umap in zip(genelabels,r_gene_umap,p_gene_umap)]

# pull out all genes which have met r and p value thresholds
def max_cor_genes(gene_list, r_thresh, p_thresh):
    if r_thresh >= 0:
        return [gene_list[i]['name'] for i in range(len(gene_list)) 
               if gene_list[i]['r_gene'] >= r_thresh and gene_list[i]['p_gene'] < p_thresh]
    if r_thresh < 0:
        return [gene_list[i]['name'] for i in range(len(gene_list)) 
               if gene_list[i]['r_gene'] <= r_thresh and gene_list[i]['p_gene'] < p_thresh]


gene_pca_pos = max_cor_genes(gene_pca_list,0.5,0.05)
gene_pca_neg = max_cor_genes(gene_pca_list,-0.5,0.05)

gene_umap_pos = max_cor_genes(gene_umap_list,0.5,0.05)
gene_umap_neg = max_cor_genes(gene_umap_list,-0.5,0.05)

# Write list into a text file

def store_list(gene_list, name):
    with open(name + str('.txt'), 'w') as f:
        for line in gene_list:
            f.write("%s\n" % line)

neg_combined = list(set(gene_umap_neg).intersection(gene_pca_neg))

pos_combined = list(set(gene_umap_pos).intersection(gene_pca_pos))

store_list(neg_combined,"significant_negative")

store_list(pos_combined,"significant_positive")