from cross_disorder_copy import cross_disorder_effect_z
import umap.plot
import matplotlib.pyplot as plt


# UMAP
umap_cor, umap_sub, umap_comp, names = cross_disorder_effect_z(method="umap")


# Cortex UMAP plot
corplot = umap.plot.points(umap_cor)
corplot.set_title('Cortex')
plt.show()

# Subcortex UMAP plot
subplot = umap.plot.points(umap_sub)
subplot.set_title('Subcortex')
plt.show()

components, variance, names = cross_disorder_effect_z(method='pca')



plt.scatter(components['cortex'][:,0], components['cortex'][:,1],color="#ff7f0e")
plt.title('Cortex Component 1 v Component 2')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.show()

plt.scatter(components['subcortex'][:,0], components['subcortex'][:,1])
plt.title('Subcortex Component 1 v Component 2')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.show()

# PCA first component vs UMAP first component
plt.scatter(components['cortex'][:,1], umap_comp['cortex'][:,1])
plt.title('PCA 2nd Components vs UMAP 2nd component')
plt.xlabel('PCA')
plt.ylabel('UMAP')
plt.show()
