import numpy as np
from sklearn.decomposition import PCA
from scipy.stats import zscore

import umap
import umap.plot

from enigmatoolbox.datasets.base import load_summary_stats


def cross_disorder_effect_z(disorder='all_disorder', measure=None,
                            additional_data_cortex=None, additional_name_cortex=None, additional_data_subcortex=None,
                            additional_name_subcortex=None, ignore=None, include=None, method='pca'):
    """Cross-disorder effect (authors: @boyongpark, @saratheriver)

        Parameters
        ----------
        disorder : list, optional
            Any combination of disorder name. Default is all available disorders, except 'adhd'. Options are:
            {'22q', 'adhd', 'asd', 'bipolar', 'depression', 'epilepsy', 'ocd', 'schizophrenia'}.
        measure : list, optional
            Any combination of measure names. Default is {'CortThick', 'CortSurf', 'SubVol'}.
        additional_data_cortex : ndarray, optional
            Name for additional cortical ENIGMA-type data. Must also provide 'additional_name_cortex'.
        additional_name_cortex : list, optional
            Additional cortical ENIGMA-type data (n, 68). Must also provide 'additional_name_cortex'.
        additional_data_subcortex : ndarray, optional
            Name for additional subcortical ENIGMA-type data. Must also provide 'additional_name_subcortex'.
        additional_name_subcortex : list, optional
            Additional subcortical ENIGMA-type data (n, 16). Must also provide 'additional_name_subcortex'.
        ignore : list, optional
            Ignore summary statistics with these expressions. Default is ('mega') as it contains NaNs.
        include : list, optional
            Include only summary statistics with these expressions. Default is empty, i.e., include everything.
        method : string, optional
            Analysis method {'pca', 'correlation','umap'}. Default is 'pca'.

        Returns
        -------
        components : dict
            Principal components of shared effects in descending order in terms of component variance.
            Only if method is 'pca'.
        variance : dict
            Variance of components. Only if method is 'pca'.
        map_cor: UMAP fit of cortex data
        , only if method is 'umap'
        map_sub: UMAP fit of subcortex data, only if method is 'umap'
        correlation_matrix : dict
            Correlation matrices of for every pair of shared effect maps. Only is method is 'correlation'.
        names : dict
            Names of disorder and case-control effect maps included in analysis.
    """
    if measure is None:
        measure = ['CortThick', 'CortSurf', 'SubVol']
    if ignore is None:
        ignore = ['mega']
    if include is None:
        include = []
    if disorder is 'all_disorder':
        disorder = ['22q', 'asd', 'bipolar', 'depression', 'epilepsy', 'ocd', 'schizophrenia']

    mat_d = {'cortex': [], 'subcortex': []}
    names = {'cortex': [], 'subcortex': []}
    for _, ii in enumerate(disorder):
        # Load summary statistics
        sum_stats = load_summary_stats(ii)
        fieldos = list(sum_stats.keys())

        # Loop through structure fields (case-control options)
        for _, jj in enumerate(fieldos):
            if 'Cort' in jj:
                if not include:
                    if not any(ig in jj for ig in ignore) and any(meas in jj for meas in measure):
                        mat_d['cortex'].append(sum_stats[jj].iloc[:, 2])
                        names['cortex'].append(ii + ': ' + jj)

                elif include:
                    if any(inc in jj for inc in include) and not any(ig in jj for ig in ignore) \
                            and any(meas in jj for meas in measure):
                        mat_d['cortex'].append(sum_stats[jj].iloc[:, 2])
                        names['cortex'].append(ii + ': ' + jj)

            if 'Sub' in jj:
                if not include:
                    if not any(ig in jj for ig in ignore) and any(meas in jj for meas in measure):
                        mat_d['subcortex'].append(sum_stats[jj].iloc[:, 2])
                        names['subcortex'].append(ii + ': ' + jj)

                elif include:
                    if any(inc in jj for inc in include) and not any(ig in jj for ig in ignore) \
                            and any(meas in jj for meas in measure):
                        mat_d['subcortex'].append(sum_stats[jj].iloc[:, 2])
                        names['subcortex'].append(ii + ': ' + jj)

    for ii, jj in enumerate(mat_d):
        mat_d[jj] = (np.asarray(mat_d[jj]))

    # If additional data and name
    if additional_data_cortex is not None and additional_name_cortex is not None:
        mat_d['cortex'] = np.append(mat_d['cortex'], additional_data_cortex)
        names['cortex'] = np.append(names['cortex'], additional_name_cortex)

    if additional_data_subcortex is not None and additional_name_subcortex is not None:
        mat_d['subcortex'] = np.append(mat_d['subcortex'], additional_data_subcortex)
        names['subcortex'] = np.append(names['subcortex'], additional_name_subcortex)

    # Z score and transpose before applying dimension reduction
    new = list()
    for i in (mat_d['cortex']):
        x = zscore(i)
        new.append(x.copy())
    new = np.asarray(new)
    mat_d['cortex']=new
    mat_d['subcortex']=zscore(mat_d['subcortex'])

    mat_d['cortex'] = np.transpose(mat_d['cortex'])
    mat_d['subcortex'] = np.transpose(mat_d['subcortex'])

    mat_d['cortex'][np.isnan(mat_d['cortex'])] = 0
    mat_d['subcortex'][np.isnan(mat_d['subcortex'])] = 0

    if method == 'pca':
        components = {'cortex': [], 'subcortex': []}
        variance = {'cortex': [], 'subcortex': []}

        pca = PCA()
        components['cortex'] = pca.fit_transform(mat_d['cortex'])
        variance['cortex'] = pca.explained_variance_ratio_

        components['subcortex'] = pca.fit_transform(mat_d['subcortex'])
        variance['subcortex'] = pca.explained_variance_ratio_

        return components, variance, names

    elif method == 'correlation':
        # noinspection PyDictCreation
        correlation_matrix = {'cortex': [], 'subcortex': []}
        correlation_matrix['cortex'] = np.corrcoef(np.transpose(mat_d['cortex'])) #transpose back to original axes
        correlation_matrix['subcortex'] = np.corrcoef(np.transpose(mat_d['subcortex']))

        return correlation_matrix, names

    elif method == "umap":
        mapper_cor = umap.UMAP(n_components=2,random_state=42).fit((mat_d['cortex']))
        mapper_sub = umap.UMAP(n_components=2,random_state=42).fit((mat_d['subcortex']))
        umap_comp = {'cortex': mapper_cor.embedding_,
                     'subcortex':mapper_sub.embedding_}
        return mapper_cor, mapper_sub, umap_comp, names

    elif method == "umap3":
        mapper_cor = umap.UMAP(n_components=3,random_state=42).fit((mat_d['cortex']))
        mapper_sub = umap.UMAP(n_components=3,random_state=42).fit((mat_d['subcortex']))
        umap_comp = {'cortex': mapper_cor.embedding_,
                     'subcortex':mapper_sub.embedding_}
        return mapper_cor, mapper_sub, umap_comp, names 