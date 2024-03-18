from enigmatoolbox.cross_disorder import cross_disorder_effect
from numpy import corrcoef
from seaborn import heatmap
import matplotlib.pyplot as plt

from cross_disorder_copy import cross_disorder_effect_z
import enigmatoolbox.cross_disorder

# Extract shared disorder effects
components, variance, names = enigmatoolbox.cross_disorder.cross_disorder_effect()
components_z, variance_z, names_z = cross_disorder_effect_z()

print(variance['cortex'][:5])

# Plot correlation matrices between zscored and default PCA output components (10 largest components
heatmap(corrcoef(components['cortex'][:10], components_z['cortex'][:10]))

plt.title('Correlation matrix of top 10 PCA components in Z-scored vs default data', fontsize = 20)
plt.xlabel('Default', fontsize = 15)
plt.ylabel('Z-scored', fontsize = 15)
plt.show()

# Find loading matrix by multiplying each component by the sqare root of its corresponding eigenvalue
