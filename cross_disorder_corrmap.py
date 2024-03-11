from enigmatoolbox.cross_disorder import cross_disorder_effect
from numpy import corrcoef
from seaborn import heatmap

# Extract shared disorder effects

correlation_matrix, names = cross_disorder_effect(method='correlation')

# Plot correlation matrices



heatmap(corrcoef(components['cortex'], components_z['cortex']), annot=True)