from enigmatoolbox.plotting import plot_subcortical, plot_cortical
import matplotlib.pyplot as plt
from enigmatoolbox.utils import parcel_to_surface
from cross_disorder_copy import cross_disorder_effect_z
import numpy as np

# Visualize the first cortical and subcortical components on the surface brains

components_z, variance_z, names_z = cross_disorder_effect_z(measure='CortThick')
umap_cor, umap_sub, umap_comp, names = cross_disorder_effect_z(measure='CortThick', method="umap")

print(umap_comp['cortex'].shape)

print(components_z['cortex'].shape)

# Z scored PCA brain maps
plot_cortical(parcel_to_surface(components_z['cortex'][:, 0], 'aparc_fsa5'), color_range=(-0.5, 0.5),

              cmap='RdBu_r', color_bar=True, size=(800, 400))

#plt.show()

plot_subcortical(components_z['subcortex'][:, 0], color_range=(-0.5, 0.5),

                 cmap='RdBu_r', color_bar=True, size=(800, 400))

#plt.show()

# Z scored UMAP brain maps
plot_cortical(parcel_to_surface(umap_comp['cortex'][:, 0], 'aparc_fsa5'), color_range=(-0.5, 0.5),

              cmap='RdBu_r', color_bar=True, size=(800, 400))

plt.show()

plot_subcortical(umap_comp['subcortex'][:, 0], color_range=(-0.5, 0.5),

                 cmap='RdBu_r', color_bar=True, size=(800, 400))

plt.show()