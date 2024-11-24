"""
Venn Diagram codes (using matplotlib-venn)
"""
from matplotlib_venn import venn2
from matplotlib.pyplot import show,title,xticks
import numpy as np
import pandas as pd
from seaborn import barplot

umap_spin = np.loadtxt('significant.txt',dtype='str')
pca_spin = np.loadtxt('significantpcs.txt',dtype='str')
umap = 5762
pca = 5601
both = len(set.intersection(set(umap_spin),set(pca_spin)))


venn2(subsets=(umap-both,pca-both,both),set_labels=("UMAP","PCA"))

venn2(subsets=(umap-both,pca-both,both),set_labels=("UMAP","PCA"))
title("Surviving Genes after Reduction")
show()


umap_up = np.loadtxt("gtex_v8_ts_DEGumap.txt",usecols=(1,3,5),skiprows=1,dtype={'names':('Tissue','# of genes',"adjusted P value"),'formats':('S50','i4',"float64")})[:54]
pca_up = np.loadtxt("gtex_v8_ts_DEGpca.txt",usecols=(1,3,5),skiprows=1,dtype={'names':('Tissue','# of genes',"adjusted P value"),'formats':('S50','i4',"float64")})[:54]

umap_down = np.loadtxt("gtex_v8_ts_DEGumap.txt",usecols=(1,3,5),skiprows=55,dtype={'names':('Tissue','# of genes',"adjusted P value"),'formats':('S50','i4',"float64")})[:54]
pca_down = np.loadtxt("gtex_v8_ts_DEGpca.txt",usecols=(1,3,5),skiprows=55,dtype={'names':('Tissue','# of genes',"adjusted P value"),'formats':('S50','i4',"float64")})[:54]

umap_both = np.loadtxt("gtex_v8_ts_DEGumap.txt",usecols=(1,3,5),skiprows=109,dtype={'names':('Tissue','# of genes',"adjusted P value"),'formats':('S50','i4',"float64")})[:54]
pca_both = np.loadtxt("gtex_v8_ts_DEGpca.txt",usecols=(1,3,5),skiprows=109,dtype={'names':('Tissue','# of genes',"adjusted P value"),'formats':('S50','i4',"float64")})[:54]


labels = []
up_umap = []
up_pca = []
d_umap =[]
d_pca = []
b_umap = []
b_pca = []

for i in umap_up:
    labels.append(i[0])
    up_umap.append(i[1])

for i in pca_up:
    up_pca.append(i[1])

for i in pca_down:
    d_pca.append(i[1])

for i in umap_down:
    d_umap.append(i[1])

for i in pca_both:
    b_pca.append(i[1])

for i in umap_both:
    b_umap.append(i[1])


dfpca=pd.DataFrame({'x':labels,'y':up_pca})
dfumap=pd.DataFrame({'x':labels,'y':up_umap})
dfpca['Reduction']="PCA"
dfumap['Reduction']="UMAP"
res=pd.concat([dfpca,dfumap])
barplot(x='x',y='y',data=res,hue='Reduction')
title("Tissue Specifity of Upregulated Genes")
xticks(rotation=90)
show()

dfpca=pd.DataFrame({'x':labels,'y':d_pca})
dfumap=pd.DataFrame({'x':labels,'y':d_umap})
dfpca['Reduction']="PCA"
dfumap['Reduction']="UMAP"
res=pd.concat([dfpca,dfumap])
barplot(x='x',y='y',data=res,hue='Reduction')
title("Tissue Specifity of Downregulated Genes")
xticks(rotation=90)
show()

dfpca=pd.DataFrame({'x':labels,'y':b_pca})
dfumap=pd.DataFrame({'x':labels,'y':b_umap})
dfpca['Reduction']="PCA"
dfumap['Reduction']="UMAP"
res=pd.concat([dfpca,dfumap])
barplot(x='x',y='y',data=res,hue='Reduction')
title("Tissue Specifity of Two-sided DEGs")
xticks(rotation=90)
show()




umap_up =['Adipose_Subcutaneous',
'Adipose_Visceral_Omentum',
'Adrenal_Gland',
'Artery_Aorta',
'Artery_Coronary',
'Artery_Tibial',
'Bladder',
'Brain_Amygdala',
'Brain_Anterior_cingulate_cortex_BA24',
'Brain_Caudate_basal_ganglia',
'Brain_Cerebellar_Hemisphere',
'Brain_Cerebellum',
'Brain_Cortex',
'Brain_Frontal_Cortex_BA9',
'Brain_Hippocampus',
'Brain_Hypothalamus',
'Brain_Nucleus_accumbens_basal_ganglia',
'Brain_Putamen_basal_ganglia',
'Brain_Spinal_cord_cervical_c-1',
'Brain_Substantia_nigra',
'Breast_Mammary_Tissue',
'Cells_Cultured_fibroblasts',
'Cells_EBVtransformed_lymphocytes',
'Cervix_Ectocervix',
'Cervix_Endocervix',
'Colon_Sigmoid',
'Esophagus_Gastroesophageal_Junction',
'Esophagus_Mucosa',
'Esophagus_Muscularis',
'Fallopian_Tube',
'Heart_Atrial_Appendage',
'Heart_Left_Ventricle',
'Kidney_Medulla',
'Liver',
'Lung',
'Muscle_Skeletal',
'Nerve_Tibial',
'Ovary',
'Pituitary',
'Skin_Not_Sun_Exposed_Suprapubic',
'Skin_Sun_Exposed_Lower_leg',
'Thyroid',
'Uterus',
'Vagina']
umap_down=['Adipose_Subcutaneous',
'Adipose_Visceral_Omentum',
'Adrenal_Gland',
'Artery_Aorta',
'Artery_Coronary',
'Artery_Tibial',
'Bladder',
'Brain_Amygdala',
'Brain_Anterior_cingulate_cortex_BA24',
'Brain_Caudate_basal_ganglia',
'Brain_Cerebellar_Hemisphere',
'Brain_Cerebellum',
'Brain_Cortex',
'Brain_Frontal_Cortex_BA9',
'Brain_Hippocampus',
'Brain_Hypothalamus',
'Brain_Nucleus_accumbens_basal_ganglia',
'Brain_Putamen_basal_ganglia',
'Brain_Spinal_cord_cervical_c-1',
'Brain_Substantia_nigra',
'Breast_Mammary_Tissue',
'Cells_Cultured_fibroblasts',
'Cells_EBVtransformed_lymphocytes',
'Colon_Sigmoid',
'Colon_Transverse',
'Esophagus_Gastroesophageal_Junction',
'Esophagus_Mucosa',
'Esophagus_Muscularis',
'Heart_Atrial_Appendage',
'Heart_Left_Ventricle',
'Kidney_Cortex',
'Liver',
'Lung',
'Minor_Salivary_Gland',
'Muscle_Skeletal',
'Nerve_Tibial',
'Ovary',
'Pancreas',
'Pituitary',
'Prostate',
'Skin_Not_Sun_Exposed_Suprapubic',
'Skin_Sun_Exposed_Lower_leg',
'Small_Intestine_Terminal_Ileum',
'Spleen',
'Stomach',
'Testis',
'Thyroid',
'Uterus',
'Vagina',
'Whole_Blood']
umap_both=['Adipose_Subcutaneous',
'Adipose_Visceral_Omentum',
'Adrenal_Gland',
'Artery_Aorta',
'Artery_Coronary',
'Artery_Tibial',
'Bladder',
'Brain_Amygdala',
'Brain_Anterior_cingulate_cortex_BA24',
'Brain_Caudate_basal_ganglia',
'Brain_Cerebellar_Hemisphere',
'Brain_Cerebellum',
'Brain_Cortex',
'Brain_Frontal_Cortex_BA9',
'Brain_Hippocampus',
'Brain_Hypothalamus',
'Brain_Nucleus_accumbens_basal_ganglia',
'Brain_Putamen_basal_ganglia',
'Brain_Spinal_cord_cervical_c-1',
'Brain_Substantia_nigra',
'Breast_Mammary_Tissue',
'Cells_Cultured_fibroblasts',
'Cells_EBVtransformed_lymphocytes',
'Cervix_Ectocervix',
'Cervix_Endocervix',
'Colon_Sigmoid',
'Colon_Transverse',
'Esophagus_Gastroesophageal_Junction',
'Esophagus_Mucosa',
'Esophagus_Muscularis',
'Fallopian_Tube',
'Heart_Atrial_Appendage',
'Heart_Left_Ventricle',
'Kidney_Cortex',
'Kidney_Medulla',
'Liver',
'Lung',
'Minor_Salivary_Gland',
'Muscle_Skeletal',
'Nerve_Tibial',
'Ovary',
'Pancreas',
'Pituitary',
'Prostate',
'Skin_Not_Sun_Exposed_Suprapubic',
'Skin_Sun_Exposed_Lower_leg',
'Small_Intestine_Terminal_Ileum',
'Spleen',
'Stomach',
'Thyroid',
'Uterus',
'Vagina',
'Whole_Blood']

pca_up=['Adipose_Subcutaneous',
'Adrenal_Gland',
'Artery_Aorta',
'Artery_Coronary',
'Artery_Tibial',
'Brain_Amygdala',
'Brain_Anterior_cingulate_cortex_BA24',
'Brain_Caudate_basal_ganglia',
'Brain_Cerebellar_Hemisphere',
'Brain_Cerebellum',
'Brain_Cortex',
'Brain_Frontal_Cortex_BA9',
'Brain_Hippocampus',
'Brain_Hypothalamus',
'Brain_Nucleus_accumbens_basal_ganglia',
'Brain_Putamen_basal_ganglia',
'Brain_Spinal_cord_cervical_c-1',
'Brain_Substantia_nigra',
'Breast_Mammary_Tissue',
'Cells_Cultured_fibroblasts',
'Cells_EBVtransformed_lymphocytes',
'Cervix_Endocervix',
'Colon_Sigmoid',
'Fallopian_Tube',
'Kidney_Medulla',
'Nerve_Tibial',
'Ovary',
'Pituitary',
'Small_Intestine_Terminal_Ileum',
'Uterus',
'Thyroid']

pca_down=['Adipose_Subcutaneous',
'Adipose_Visceral_Omentum',
'Adrenal_Gland',
'Artery_Aorta',
'Artery_Coronary',
'Artery_Tibial',
'Brain_Amygdala',
'Brain_Anterior_cingulate_cortex_BA24',
'Brain_Caudate_basal_ganglia',
'Brain_Cerebellar_Hemisphere',
'Brain_Cerebellum',
'Brain_Cortex',
'Brain_Frontal_Cortex_BA9',
'Brain_Hippocampus',
'Brain_Hypothalamus',
'Brain_Nucleus_accumbens_basal_ganglia',
'Brain_Putamen_basal_ganglia',
'Brain_Spinal_cord_cervical_c-1',
'Brain_Substantia_nigra',
'Breast_Mammary_Tissue',
'Cells_Cultured_fibroblasts',
'Cells_EBVtransformed_lymphocytes',
'Colon_Sigmoid',
'Colon_Transverse',
'Esophagus_Gastroesophageal_Junction',
'Esophagus_Mucosa',
'Esophagus_Muscularis',
'Heart_Atrial_Appendage',
'Heart_Left_Ventricle',
'Kidney_Cortex',
'Liver',
'Lung',
'Minor_Salivary_Gland',
'Muscle_Skeletal',
'Nerve_Tibial',
'Ovary',
'Pancreas',
'Pituitary',
'Prostate',
'Skin_Not_Sun_Exposed_Suprapubic',
'Skin_Sun_Exposed_Lower_leg',
'Small_Intestine_Terminal_Ileum',
'Spleen',
'Stomach',
'Testis',
'Thyroid',
'Uterus',
'Vagina',
'Whole_Blood']


pca_both=['Adipose_Subcutaneous',
'Adipose_Visceral_Omentum',
'Adrenal_Gland',
'Artery_Aorta',
'Artery_Coronary',
'Artery_Tibial',
'Bladder',
'Brain_Amygdala',
'Brain_Anterior_cingulate_cortex_BA24',
'Brain_Caudate_basal_ganglia',
'Brain_Cerebellar_Hemisphere',
'Brain_Cerebellum',
'Brain_Cortex',
'Brain_Frontal_Cortex_BA9',
'Brain_Hippocampus',
'Brain_Hypothalamus',
'Brain_Nucleus_accumbens_basal_ganglia',
'Brain_Putamen_basal_ganglia',
'Brain_Spinal_cord_cervical_c-1',
'Brain_Substantia_nigra',
'Breast_Mammary_Tissue',
'Cells_Cultured_fibroblasts',
'Cells_EBVtransformed_lymphocytes',
'Cervix_Endocervix',
'Colon_Sigmoid',
'Colon_Transverse',
'Esophagus_Gastroesophageal_Junction',
'Esophagus_Mucosa',
'Esophagus_Muscularis',
'Fallopian_Tube',
'Heart_Atrial_Appendage',
'Heart_Left_Ventricle',
'Kidney_Cortex',
'Liver',
'Minor_Salivary_Gland',
'Muscle_Skeletal',
'Nerve_Tibial',
'Ovary',
'Pancreas',
'Pituitary',
'Skin_Not_Sun_Exposed_Suprapubic',
'Skin_Sun_Exposed_Lower_leg',
'Spleen',
'Stomach',
'Thyroid',
'Uterus',
'Whole_Blood']


'''
venn2((set(umap_up), set(pca_up)),set_labels=("UMAP","PCA"))
title("Tissue Specificity of Upregulated Genes")
show()

venn2((set(umap_down), set(pca_down)),set_labels=("UMAP","PCA"))
title("Tissue Specificity of Downregulated Genes")
show()

venn2((set(umap_both), set(pca_both)),set_labels=("UMAP","PCA"))
title("Tissue Specificity of Regulated Genes")
show()
'''
