# gene expression analysis
from enigmatoolbox.datasets import fetch_ahba
import numpy as np
from enigmatoolbox.utils.parcellation import parcel_to_surface, surface_to_parcel
from scipy.stats import spearmanr, false_discovery_control
from enigmatoolbox.permutation_testing import spin_test
from statsmodels.stats.multitest import multipletests

# correlate PC1 with gene expression data from brain regions using spin test
genes = fetch_ahba()
reglabels = genes['label']
genelabels = list(genes.columns)

genes = fetch_ahba()
reglabels = genes['label']

# remove the label column
genes = genes.drop('label', axis=1)
genelabels = list(genes.columns)

def corr_gene_exp(solution):
    if len(solution) == 200:
        print("Schaefer")
        # put it on a surface
        d_mri_surf = parcel_to_surface(solution, f'schaefer_200_conte69')
        # map it back on DK parcels used by ENIGMA
        d_mri = surface_to_parcel(d_mri_surf, f'aparc_conte69')
        # remove midbrain
        d_mri = np.delete(d_mri, np.array([0, 4, 39]))
    elif len(solution) == 68:
        print("DK")
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

def corr_gene_exp_fdr(solution):
    if len(solution) == 200:
        # put it on a surface
        d_mri_surf = parcel_to_surface(solution, f'schaefer_200_conte69')
        # map it back on DK parcels used by ENIGMA
        d_mri = surface_to_parcel(d_mri_surf, f'aparc_conte69')
        # remove midbrain
        d_mri = np.delete(d_mri, np.array([0, 4, 39]))
    elif len(solution) == 68:
        d_mri = solution #parcellations as denoted by the first or second component of PCA or UMAP
    else:
        raise ValueError('Solution must be 200 or 68 long')
    r_gene = [] #R is correlation coefficient
    p_gene = [] #p is probability that correlation exists if null is true
    for i in enumerate(genelabels):
        geneexp = genes.iloc[0:68,i]
        # get the location of non-nan values
        idx = ~np.isnan(geneexp)   # idx is a list of locations    
        c_gene = spearmanr(d_mri[idx], geneexp[idx])
        r_gene = np.append(r_gene, c_gene[0])
        p_gene = np.append(p_gene, c_gene[1]) 
    # adjust p values with FDR (Benjami-Hochberg)
    p_gene = false_discovery_control(p_gene)
    
    return r_gene, p_gene


def corr_gene_exp_spin(solution):
    d_mri = solution #parcellations as denoted by the first or second component of PCA or UMAP
    p_values = np.zeros(len(genelabels))
    for i in range(len(genelabels)):
        geneexp = genes[genelabels[i]]
        geneexp = geneexp[0:68]
        p_values[i] = spin_test(geneexp, d_mri, surface_name='fsa5', parcellation_name='aparc',
                        type='spearman', n_rot=1000, null_dist=False)  
    # apply Benjamini-Hochberg correction

    p_values_corrected = multipletests(p_values, alpha=0.05, method='fdr_bh')[1]
    # print how many genes are significant
    return(f'Number of significant genes: {np.sum(p_values_corrected < 0.05)}')
    
# pull out all genes which have met r and p value thresholds
def max_cor_genes(gene_list, p_thresh):
    return [gene_list[i]['name'] for i in range(len(gene_list)) 
        if gene_list[i]['p_gene'] < p_thresh]




# Write list into a text file

def store_list(gene_list, name):
    with open(name + str('.txt'), 'w') as f:
        for line in gene_list:
            f.write("%s\n" % line)
