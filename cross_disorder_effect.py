from cross_disorder_copy import cross_disorder_effect

from enigmatoolbox.plotting import plot_cortical, plot_subcortical

from enigmatoolbox.utils import parcel_to_surface

import matplotlib.pyplot as plt

# Extract shared disorder effects

components, variance, names = cross_disorder_effect()

# Visualize cortical and subcortical eigenvalues in scree plots

fig, ax = plt.subplots(1, 2, figsize=(14, 6))

for ii, jj in enumerate(components):

   ax[ii].plot(variance[jj], lw=2, color='#A8221C', zorder=1)

   ax[ii].scatter(range(variance[jj].size), variance[jj], s=78, color='#A8221C',

                  linewidth=1.5, edgecolor='w', zorder=3)

   ax[ii].set_xlabel('Components')

   if ii == 0:

        ax[ii].set_ylabel('Cortical eigenvalues')

   else:

        ax[ii].set_ylabel('Subcortical eigenvalues')

   ax[ii].spines['top'].set_visible(False)

   ax[ii].spines['right'].set_visible(False)

fig.tight_layout()

plt.show()

# Visualize the first cortical and subcortical components on the surface brains

plot_cortical(parcel_to_surface(components['cortex'][:, 0], 'aparc_fsa5'), color_range=(-0.5, 0.5),

              cmap='RdBu_r', color_bar=True, size=(800, 400))

plot_subcortical(components['subcortex'][:, 0], color_range=(-0.5, 0.5),

                 cmap='RdBu_r', color_bar=True, size=(800, 400))