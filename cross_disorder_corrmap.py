from enigmatoolbox.cross_disorder import cross_disorder_effect
from numpy import corrcoef, sqrt, concatenate, newaxis
from seaborn import heatmap
import matplotlib.pyplot as plt


from cross_disorder_copy import cross_disorder_effect_z
import enigmatoolbox.cross_disorder

# Extract shared disorder effect PCA
components_z, variance_z, names_z = cross_disorder_effect_z(measure='CortThick')

#UMAP
umap_cor, umap_sub, umap_comp, names = cross_disorder_effect_z(measure='CortThick',method="umap")

heatmap(corrcoef(components_z['cortex'], umap_comp['cortex'][0]))
plt.title('PCA components vs UMAP components')
plt.xlabel('PCA components')
plt.ylabel('UMAP components')
plt.show()

# Get correlation matrix
corr_matrix, names_m = enigmatoolbox.cross_disorder.cross_disorder_effect(method="correlation")

# Plot correlation matrices between zscored and original variables (10 largest components)
#heatmap(corrcoef(components_z['cortex'][:10], corr_matrix['cortex']), yticklabels=names['cortex'])

# plt.title('Cortex Correlation matrix vs Original Variables', fontsize=20)
#plt.ylabel('Original Variables', fontsize=15)
#plt.xlabel('PCA components from Z-Scored', fontsize=15)
#plt.show()

# Plot correlation matrices between zscored and original variables (10 largest components)
#heatmap(corrcoef(components['cortex'][:10], corr_matrix['cortex']), yticklabels=names['cortex'])

# plt.title('Cortex Correlation matrix vs Original Variables', fontsize=20)
#plt.ylabel('Original Variables', fontsize=15)
#plt.xlabel('PCA components', fontsize=15)
#plt.show()

# Find loading matrix by multiplying each component by the square root of its corresponding eigenvalue
#load = components * sqrt(variance['cortex'])

# load for z scored

# load_z = components_z * sqrt(variance_z['cortex'])

