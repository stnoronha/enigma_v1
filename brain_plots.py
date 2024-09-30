from enigmatoolbox.plotting import plot_subcortical, plot_cortical
import matplotlib.pyplot as plt
from enigmatoolbox.utils import parcel_to_surface
from cross_disorder_copy import cross_disorder_effect_z
import numpy as np
from enigmatoolbox.datasets.base import load_summary_stats
from enigmatoolbox.utils.parcellation import parcel_to_surface


# Visualize the first cortical and subcortical components on the surface brains

components_z, variance_z, names_z = cross_disorder_effect_z(measure="CortThick")
umap_cor, umap_sub, umap_comp, names = cross_disorder_effect_z(measure="CortThick",method="umap")

print(names)
print(components_z['cortex'][:, 0])

# Z scored PCA brain maps for z-scored first component
plot_cortical(parcel_to_surface(components_z['cortex'][:, 0], 'aparc_fsa5'), color_range=(-5, 6),
              cmap='RdBu_r', color_bar=True, size=(800, 400))

plt.show()


# Z scored UMAP brain maps
plot_cortical(parcel_to_surface(umap_comp['cortex'][:, 0], 'aparc_fsa5'), surface_name='fsa5', color_range=(-2, 8),
              cmap='RdBu_r', color_bar=True, size=(800, 400))

plt.show()





