from enigmatoolbox.plotting import plot_subcortical, plot_cortical
import matplotlib.pyplot as plt
from enigmatoolbox.utils import parcel_to_surface
from cross_disorder_copy import cross_disorder_effect_z
# Visualize the first cortical and subcortical components on the surface brains

components_z, variance_z, names_z = cross_disorder_effect_z()

# Z scored PCA brain maps
plot_cortical(parcel_to_surface(components_z['cortex'][:, 0], 'aparc_fsa5'), color_range=(-0.5, 0.5),

              cmap='RdBu_r', color_bar=True, size=(800, 400))

plt.show()

plot_subcortical(components_z['subcortex'][:, 0], color_range=(-0.5, 0.5),

                 cmap='RdBu_r', color_bar=True, size=(800, 400))

plt.show()
