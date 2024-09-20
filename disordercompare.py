from enigmatoolbox.cross_disorder import cross_disorder_effect
from cross_disorder_copy import cross_disorder_effect_z
from numpy import array,transpose
from seaborn import heatmap, barplot
import matplotlib.pyplot as plt
import pandas as pd



disorder_list= ['22q', 'asd', 'bipolar', 'depression', 'epilepsy', 'ocd', 'schizophrenia']
explained = {}
for d in disorder_list:
    # PCA analysis of cross_disorder for cortical thickness
    components, variance, names=cross_disorder_effect(disorder=[d],measure="CortThick")
    # get explained variance ratio for 
    v=variance['cortex'][0]
    explained[d]=v
expdf = pd.DataFrame(list(explained.keys()),list(explained.values()))
barplot(expdf)
plt.show()


explained_z = {}
for d in disorder_list:
    # PCA analysis of cross_disorder for cortical thickness
    components, variance, names=cross_disorder_effect_z(disorder=[d],measure="CortThick")
    # get explained variance ratio for 
    v = variance['cortex'][0]
    #if len(v) == 1:
    #    print(v)
    #    explained_z[d]=v
    explained_z[d]=v
barplot(explained_z)
heatmap(explained_z,annot=True,square=True,yticklabels=["Component 1","Component 2"],xticklabels=disorder_list)
plt.title("Disorder Explained Variance Ratio for 1st and 2nd Z scored PCA components")
plt.show()