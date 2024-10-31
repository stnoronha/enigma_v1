
from cross_disorder_copy import cross_disorder_effect_z
from enigmatoolbox.cross_disorder import cross_disorder_effect
from gene_expression import store_list
from enigmatoolbox.utils.parcellation import parcel_to_surface
from enigmatoolbox.plotting import plot_cortical
from statistics import median

# Visualize the first cortical and components on the surface brains

components_z, variance_z, names_z = cross_disorder_effect_z(measure="CortThick")
components,variance,names = cross_disorder_effect(measure="CortThick")
umap_cor, umap_sub, umap_comp, names = cross_disorder_effect_z(measure="CortThick",method="umap3")

store_list(components['cortex'][:, 0],"PCAnon-dk")
store_list(components_z['cortex'][:, 0],"PCA-dk")
store_list(umap_comp['cortex'][:, 0],"UMAP3-dk")

med_pca1 = median(components['cortex'][:, 0])
med_pca = median(components_z['cortex'][:, 0])
med_umap = median(umap_comp['cortex'][:, 0])

range_pca1 = (max(components['cortex'][:, 0])-min(components['cortex'][:, 0]))/2.0
range_pca = (max(components_z['cortex'][:, 0])-min(components_z['cortex'][:, 0]))/2.0
range_umap = (max(umap_comp['cortex'][:, 0])-min(umap_comp['cortex'][:, 0]))/2.0


plot_cortical(parcel_to_surface(components['cortex'][:, 0], 'aparc_fsa5'), color_range=(med_pca1-range_pca1, med_pca1+range_pca1),
             cmap='RdBu_r', color_bar=True, size=(800, 400),interactive=False)

# Z scored PCA brain maps for z-scored first component
plot_cortical(parcel_to_surface(components_z['cortex'][:, 0], 'aparc_fsa5'), color_range=(med_pca-range_pca, med_pca+range_pca),
              cmap='RdBu_r', color_bar=True, size=(800, 400),interactive=False)


# Z scored UMAP brain maps
plot_cortical(parcel_to_surface(umap_comp['cortex'][:, 0], 'aparc_fsa5'), surface_name='fsa5', color_range=(-13,18),
              cmap='RdBu_r', color_bar=True, size=(800, 400),interactive=False)






