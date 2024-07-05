# gene expression analysis
from enigmatoolbox.datasets import fetch_ahba
import matplotlib.pyplot as plt
from enigmatoolbox.utils import parcel_to_surface
from cross_disorder_copy import cross_disorder_effect_z
import numpy as np
from enigmatoolbox.utils.parcellation import parcel_to_surface, surface_to_parcel
from scipy.stats import spearmanr, false_discovery_control
from enigmatoolbox.permutation_testing import spin_test,shuf_test

# PCA cross disorder effect
components_z, variance_z, names_z = cross_disorder_effect_z()

# UMAP cross disorder effect
umap_cor, umap_sub, umap_comp, names = cross_disorder_effect_z(method="umap")

# Isomap cross disorder effect
iso_cor, iso_sub, iso_comp, names = cross_disorder_effect_z(method="isomap")

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

def corr_gene_exp_fdr(solution):
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
    # adjust p values with FDR (Benjami-Hochberg)
    p_gene = false_discovery_control(p_gene)
    return r_gene, p_gene


def corr_gene_exp_spin(solution):
# Remove subcortical values corresponding to the ventricles
# (as we don't have connectivity values for them!)
    fc_ctx_dc = np.sum(fc_ctx, axis=0)
    sc_ctx_dc = np.sum(sc_ctx, axis=0)

    # Map parcellated data to the surface
    fc_ctx_dc_fsa5 = parcel_to_surface(fc_ctx_dc, 'aparc_fsa5')

    sc_ctx_dc_fsa5 = parcel_to_surface(sc_ctx_dc, 'aparc_fsa5')
    SV_d_noVent = SV_d.drop([np.where(SV['Structure'] == 'LLatVent')[0][0],
                            np.where(SV['Structure'] == 'RLatVent')[0][0]]).reset_index(drop=True)

    # Spin permutation testing for two cortical maps
    fc_ctx_p, fc_ctx_d = spin_test(fc_ctx_dc, CT_d, surface_name='fsa5', parcellation_name='aparc',
                                type='pearson', n_rot=1000, null_dist=True)

    sc_ctx_p, sc_ctx_d = spin_test(sc_ctx_dc, CT_d, surface_name='fsa5', parcellation_name='aparc',
                                type='pearson', n_rot=1000, null_dist=True)

    # Shuf permutation testing for two subcortical maps

    fc_sctx_p, fc_sctx_d = shuf_test(fc_sctx_dc, SV_d_noVent, n_rot=1000,
                                    type='pearson', null_dist=True)

    sc_sctx_p, sc_sctx_d = shuf_test(sc_sctx_dc, SV_d_noVent, n_rot=1000,
                                    type='spearman', null_dist=True)
    # Store p-values and null distributions
    p_and_d = {'functional cortical hubs': [fc_ctx_p, fc_ctx_d], 'functional subcortical hubs': [fc_sctx_p, fc_sctx_d],
            'structural cortical hubs': [sc_ctx_p, sc_ctx_d], 'structural subcortical hubs': [sc_sctx_p, sc_sctx_d]}

# pull out all genes which have met r and p value thresholds
def max_cor_genes(gene_list, r_thresh, p_thresh):
    if r_thresh >= 0:
        return [gene_list[i]['name'] for i in range(len(gene_list)) 
               if gene_list[i]['r_gene'] >= r_thresh and gene_list[i]['p_gene'] < p_thresh]
    if r_thresh < 0:
        return [gene_list[i]['name'] for i in range(len(gene_list)) 
               if gene_list[i]['r_gene'] <= r_thresh and gene_list[i]['p_gene'] < p_thresh]


# Write list into a text file

def store_list(gene_list, name):
    with open(name + str('.txt'), 'w') as f:
        for line in gene_list:
            f.write("%s\n" % line)
