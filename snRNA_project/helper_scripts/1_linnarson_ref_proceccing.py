import scanpy as sc
import celltypist
from celltypist import models
import numpy as np
import loompy


input_folder = '/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref/input_loom/'
loom_path=f'{input_folder}linnarsson_mouse_selection-VM_STR_Ctx_final.loom'
loom_path=f'{input_folder}L5_All.loom'
output_path='/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref/anndata_from_loom/'

import pandas as pd
import anndata as ad

# Conversion of the Loom file to anndata
def convert_loom_to_anndata(loom_file_loc, ca_list=[], ra_list=[], ca_index='Description', ra_index='Gene'):
    loom_file = loompy.connect(loom_file_loc)

    attr_lists = [ca_list, ra_list]
    
    # if attr lists are empy, keep original columns/rows
    for idx, attr_list in enumerate(attr_lists):
        if len(attr_lists[idx]) == 0:
            if idx == 0: 
                attr_lists[idx] = loom_file.ca.keys()
            elif idx == 1: 
                attr_lists[idx] = loom_file.ra.keys()
    
    # select index columns for the dataframes
    attr_indexes = [ca_index, ra_index]
    for idx, index in enumerate(attr_indexes):
        if type(index) == int:
            attr_indexes[idx] = attr_lists[idx][index]
        elif type(index) == str:
            assert index in attr_lists[idx]
    print(f'The indeces for var and obs will be assigned to {attr_indexes[0]} and {attr_indexes[1]}')
    
    # create var and obs dataframes with defined columns and indexes (indices)
    ad_attr = [pd.DataFrame(), pd.DataFrame()]
    for idx, attr_list in enumerate(attr_lists):
        for attr in attr_list:
            if idx == 0: 
                ad_attr[idx][attr] = loom_file.ca[attr]
            elif idx == 1: 
                ad_attr[idx][attr] = loom_file.ra[attr]
        ad_attr[idx].index = ad_attr[idx][attr_indexes[idx]]

    adata = ad.AnnData(X = loom_file[:, :].T, var=ad_attr[1], obs=ad_attr[0])
    loom_file.close()
    return adata

# Subsampling from the anndata (created from loom and subsetted based on subset_all list described on the line 96)
def subsample_anndata(anndata, annot_column='Description', counts=[10, 1000]):
    print(f'Dataset will be downsampled to contain between {counts[0]} and {counts[1]} cells per celltype')
    anndata_subset = anndata.copy()
    cells_to_keep = []

    for x in anndata_subset.obs[annot_column].unique():
        print(x)
        all_cells = anndata_subset.obs[anndata_subset.obs[annot_column] == x]['CellID'].to_list()
        if len(all_cells) < counts[0]:
            anndata_subset = anndata_subset[anndata_subset.obs[annot_column] != x, :]
            print(f'{x} with {len(all_cells)} cells will be dropped')
        elif len(all_cells) >= counts[0] and len(all_cells) <= counts[1]:
            cells_to_keep += all_cells
            print(f'All {len(all_cells)} cells will be used')
        elif len(all_cells) > counts[1]:
            cells_that_won_the_lottery = np.random.choice(all_cells, size=counts[1], replace=False).tolist()
            print(f'{len(cells_that_won_the_lottery)} cells will be kept out of {len(all_cells)}')
            cells_to_keep += cells_that_won_the_lottery
    
    anndata_subset = anndata_subset[anndata_subset.obs['CellID'].isin(cells_to_keep), :]
    print(anndata_subset.obs[annot_column].value_counts())
    return anndata_subset

# subset_all = ["Vascular endothelial cells, arterial", "Vascular endothelial cells, venous", "Neuronal intermidate progenitor cells", "Vascular leptomeningeal cells", "Vascular smooth muscle cells, arterial", "Myelin forming oligodendrocytes (MFOL)", "Microglia", "Inhibitory neurons, midbrain", "Dopaminergic neurons, ventral midbrain (SNc, VTA)", 'Neuronal intermidate progenitor cells', 'Excitatory neurons, midbrain', 'Non-telencephalon astrocytes, protoplasmic', 'Newly formed oligodendrocytes (NFOL)','Mature oligodendrocytes', 'Dorsal midbrain Myoc-expressing astrocyte-like']

subset_all = ["Vascular endothelial cells, arterial", "Vascular endothelial cells, venous",
 "Neuronal intermidate progenitor cells", "Vascular leptomeningeal cells", 
 "Vascular smooth muscle cells, arterial",
"Myelin forming oligodendrocytes (MFOL)", "Inhibitory neurons, midbrain",
 "Dopaminergic neurons, ventral midbrain (SNc, VTA)", 'Neuronal intermidate progenitor cells', 'Excitatory neurons, midbrain',
 'Non-telencephalon astrocytes, protoplasmic',
'Mature oligodendrocytes', 'Dorsal midbrain Myoc-expressing astrocyte-like']
len(subset_all)
# subset_all = ["Excitatory neurons, midbrain", "Myelin forming oligodendrocytes (MFOL)" , "Vascular endothelial cells, arterial", "Vascular endothelial cells, venous", "Neuronal intermidate progenitor cells", "Vascular leptomeningeal cells", "Vascular smooth muscle cells, arterial", "Microglia", "Inhibitory neurons, midbrain", "Dopaminergic neurons, ventral midbrain (SNc, VTA)", 'Neuronal intermidate progenitor cells', 'Non-telencephalon astrocytes, protoplasmic','Mature oligodendrocytes', 'Dorsal midbrain Myoc-expressing astrocyte-like']
# subset_all = ["Vascular endothelial cells, arterial", "Vascular endothelial cells, venous", "Neuronal intermidate progenitor cells", "Vascular leptomeningeal cells", "Vascular smooth muscle cells, arterial", "Microglia", "Inhibitory neurons, midbrain", "Dopaminergic neurons, ventral midbrain (SNc, VTA)", 'Neuronal intermidate progenitor cells', 'Non-telencephalon astrocytes, protoplasmic', 'Dorsal midbrain Myoc-expressing astrocyte-like']

# removed excitatory neurons, midbrain and Myelin forming oligodendrocytes (MFOL) 



input_folder = '/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref/input_loom/'
output_path = '/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref/anndata_from_loom/'
#output_folder = f'{output_path}L5_{today.strftime("%d%m%y")}/'
loom_file = loompy.connect(f'{input_folder}L5_All.loom') 
# loom_file = loompy.connect(f'{input_folder}linnarsson_mouse_selection-VM_STR_Ctx_final.loom')

########################################################################################################################
loom_file = input_folder + "L5_All.loom"
# loom_file = input_folder + "linnarsson_mouse_selection-VM_STR_Ctx_final.loom"

ca_selection = ['CellID', 'Class', 'Description', 'Region']
ra_selection = ['Gene']
loom_file = loompy.connect(loom_file)

attr_lists = [ca_selection, ra_selection]
ca_index='Description'
ra_index='Gene'
# if attr lists are empy, keep original columns/rows
for idx, attr_list in enumerate(attr_lists):
    if len(attr_lists[idx]) == 0:
        if idx == 0: 
            attr_lists[idx] = loom_file.ca.keys()
        elif idx == 1: 
            attr_lists[idx] = loom_file.ra.keys()

# select index columns for the dataframes
attr_indexes = [ca_index, ra_index]
for idx, index in enumerate(attr_indexes):
    if type(index) == int:
        attr_indexes[idx] = attr_lists[idx][index]
    elif type(index) == str:
        assert index in attr_lists[idx]
print(f'The indeces for var and obs will be assigned to {attr_indexes[0]} and {attr_indexes[1]}')

# create var and obs dataframes with defined columns and indexes (indices)
ad_attr = [pd.DataFrame(), pd.DataFrame()]
for idx, attr_list in enumerate(attr_lists):
    for attr in attr_list:
        if idx == 0: 
            ad_attr[idx][attr] = loom_file.ca[attr]
        elif idx == 1: 
            ad_attr[idx][attr] = loom_file.ra[attr]
    ad_attr[idx].index = ad_attr[idx][attr_indexes[idx]]

adata = ad.AnnData(X = loom_file[:, :].T, var=ad_attr[1], obs=ad_attr[0])
loom_file.close()

adata = convert_loom_to_anndata(loom_file, ca_selection, ra_selection)
np.unique(adata.obs['Description'])
adata_selected = adata[adata.obs.Description.isin(subset_all)] # 38618 × 27998

adata_selected = subsample_anndata(adata_selected)
output_path = '/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref/anndata_from_loom/'
adata_selected.write(output_path + "/ADATA_midbrain_subset.h5ad")
l = sc.read_h5ad(output_path + "/ADATA_midbrain_subset.h5ad")
l

np.unique(adata_selected.obs['Description'])

# Save as csv files which used to construct the ref on R for SingleR
mat = pd.DataFrame(adata_selected.X)
mat.iloc[1:3,1:3]
mat.index = list(adata_selected.obs['Description'])
mat.columns = list(adata_selected.var['Gene'])
mat = mat.transpose()
mat.head()
mat_da = mat[['Dopaminergic neurons, ventral midbrain (SNc, VTA)']]        
mat_da.head()


# mat = pd.DataFrame(adata_selected.X, columns=adata_selected.var.index, index=adata_selected.obs.index)
mat.to_csv(f'{output_path}COUNT_MATRIX_ADATA_1000.csv')
adata_selected.obs['Description'].to_csv(f'{output_path}CELLS_ADATA_1000.csv')
adata_selected.obs['CellID'].to_csv(f'{output_path}CELL_IDs_ADATA_1000.csv')

adata_selected.var['Gene'].to_csv(f'{output_path}GENES_ADATA_1000.csv')
########################################################################################################################