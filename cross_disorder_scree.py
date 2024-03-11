import numpy

from cross_disorder_copy import cross_disorder_effect_z
import enigmatoolbox.cross_disorder

from pprint import pprint
from enigmatoolbox.plotting import plot_cortical, plot_subcortical

from enigmatoolbox.utils import parcel_to_surface
from seaborn import heatmap
import matplotlib.pyplot as plt

# Extract shared disorder effects for Z score
components, variance, names = enigmatoolbox.cross_disorder.cross_disorder_effect()
components_asd, variance_asd, names = enigmatoolbox.cross_disorder.cross_disorder_effect(disorder=['ocd','asd'])

components_z, variance_z, names_z = cross_disorder_effect_z()
componentsz_asd, variancez_asd, names = cross_disorder_effect_z(disorder=['ocd'])

# Print correlation matrix between original components and zscored components for cortical
heatmap(numpy.corrcoef(components['cortex'][0:5], components_z['cortex'][0:5]),annot=True)

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

fig_z, ax_z = plt.subplots(1, 2, figsize=(14, 6))

for ii, jj in enumerate(components_z):

    ax_z[ii].plot(variance_z[jj], lw=2, color='#A8221C', zorder=1)

    ax_z[ii].scatter(range(variance_z[jj].size), variance_z[jj], s=78, color='#A8221C',

                     linewidth=1.5, edgecolor='w', zorder=3)

    ax_z[ii].set_xlabel('Components')

    if ii == 0:

        ax_z[ii].set_ylabel('Cortical eigenvalues')

    else:

        ax_z[ii].set_ylabel('Subcortical eigenvalues')

    ax_z[ii].spines['top'].set_visible(False)

    ax_z[ii].spines['right'].set_visible(False)

fig.tight_layout()

plt.show()

'''
# Visualize the first cortical and subcortical components on the surface brains

plot_cortical(parcel_to_surface(components['cortex'][:, 0], 'aparc_fsa5'), color_range=(-0.5, 0.5),

              cmap='RdBu_r', color_bar=True, size=(800, 400))

plot_subcortical(components['subcortex'][:, 0], color_range=(-0.5, 0.5),

                 cmap='RdBu_r', color_bar=True, size=(800, 400))'''
