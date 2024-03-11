from enigmatoolbox.cross_disorder import cross_disorder_effect

from nilearn import plotting

# Extract shared disorder effects

correlation_matrix, names = cross_disorder_effect(method='correlation')

# Plot correlation matrices

plotting.plot_matrix(correlation_matrix['cortex'], figure=(12, 8), labels=names['cortex'], vmax=1,

                     vmin=-1, cmap='RdBu_r', auto_fit=False)

plotting.plot_matrix(correlation_matrix['subcortex'], figure=(12, 8), labels=names['subcortex'], vmax=1,

                     vmin=-1, cmap='RdBu_r', auto_fit=False)