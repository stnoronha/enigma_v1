from cross_disorder_copy import cross_disorder_effect_z
import umap.plot
import matplotlib.pyplot as plt


# UMAP
umap_cor, umap_sub, umap_comp, names = cross_disorder_effect_z(measure='CortThick',method="umap")


# Cortex UMAP plot
corplot = umap.plot.points(umap_cor)
corplot.set_title('Cortex')
#plt.show()

# Subcortex UMAP plot
subplot = umap.plot.points(umap_sub)
subplot.set_title('Subcortex')
#plt.show()


components, variance, names = cross_disorder_effect_z(measure = 'CortThick',method='pca')

plt.scatter(components['cortex'][:,0], components['cortex'][:,1])
plt.title('Cortex Component 1 v Component 2')
plt.xlabel('PC1')
plt.ylabel('PC2')
#plt.show()

plt.scatter(components['subcortex'][:,0], components['subcortex'][:,1])
plt.title('Subcortex Component 1 v Component 2')
plt.xlabel('PC1')
plt.ylabel('PC2')
#plt.show()

# PCA vs UMAP components
plt.scatter(components['cortex'][:50,0], umap_comp['cortex'][:,0])
plt.title('PCA Components vs UMAP components')
plt.xlabel('PCA')
plt.ylabel('UMAP')
plt.show()
