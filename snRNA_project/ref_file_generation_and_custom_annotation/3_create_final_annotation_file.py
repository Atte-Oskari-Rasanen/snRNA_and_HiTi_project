import pandas as pd
import numpy as np
from collections import Counter


'''
Take in the original full annotation file and merge with the further annotated subset file by going through the 
cell_i and if a match is found between the full_ann and the da_ann or glut_ann (da/glut_ann only contains the cell_ids of the subset)

'''


def merge_annotations(full_ann, sub_ann):
  full_ann_final = full_ann.copy()
  for cell_i, cell in enumerate(sub_ann.index):
      if cell in full_ann.index:
          i = list(full_ann.index).index(cell)
          full_ann_final.iloc[i,0] = sub_ann.iloc[cell_i,0]
  return(full_ann_final)

p = "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/"

p = "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/annotations_final_nonspecific_da.csv"

p = "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/final_annotations2.csv"
p = "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/annotations_final_nonspecific_da_scifi.csv"

full_ann = pd.read_csv(p)
full_ann.rename(columns = {'Unnamed: 0':'genes'}, inplace = True)
full_ann = full_ann.set_index('genes')

#p_da = "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/DA_subtypes_scifi_meta.csv"
p_da = "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/DA_subtypes_parse_meta.csv"
p_da = "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/DA_subtypes_parse_meta_2.csv"
p_da = "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/DA_subtypes_parse_meta_2_harsh.csv"


da_ann = pd.read_csv(p_da)
da_ann.rename(columns = {'Unnamed: 0':'genes'}, inplace = True)
da_ann = da_ann.set_index('genes')

full_ann_final = merge_annotations(full_ann, da_ann)


p = "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/"
full_ann_final.to_csv(p + "final_annotations_parse_with_dasubtypes_2_harsh.csv")
# full_ann_final.to_csv(p + "final_annotations_scifi_with_dasubtypes2.csv")

np.unique(full_ann_final.iloc[:,0])
Counter(list(full_ann_final.iloc[:,0]))


full_ann_final = full_ann.copy()
for cell_i, cell in enumerate(da_ann.index):
  print(cell)
  if cell in full_ann.index:
      i = list(full_ann.index).index(cell)
      full_ann_final.iloc[i,1] = da_ann.iloc[cell_i,1]

#####################################################################

p = "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/annotations_final_specific_da_scifi.csv"
p = "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/annotations_final_specific_da_scifi.csv"

full_ann_final.rename(columns = {'x':'Annotations_da_specific'}, inplace = True)
full_ann_final.to_csv

for i, cell in enumerate(full_ann_final.iloc[:,0]):
    if cell == "Substantia nigra neurons":
        full_ann_final.iloc[i,0] = "Inhibitory neurons, midbrain"
        

Counter(list(full_ann_final.iloc[:,0]))



np.unique(full_ann_final['x'])

np.unique(full_ann_final.iloc[:,1])

final_ann = da_ann.merge(da_ann, how='left').fillna("Unknown")
np.unique(final_ann['0'])

da_ann = pd.read_csv(p)

full_ann = pd.read_csv(p + "full_ann2.csv", index_col=[0])
da_ann = pd.read_csv("/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/DA_subtypes_scifi_meta.csv")
glut_ann = pd.read_csv(p + "glut_ann2.csv", index_col=[0])

final_ann = full_ann.merge(glut_ann, how='left').fillna("Unknown")
np.unique(final_ann.iloc[:,2])

############################################################################
full_final_glut_all = merge_annotations(full_ann, glut_ann)
full_final_glut_all.to_csv(p + "final_annotations2.csv")
############################################################################
Counter(list(full_final_glut_all.iloc[:,0]))


# Drop excitatory neurons annotation
full_ann_sub = full_ann[full_ann.pbmc.merged.final@meta.data$annotations_1000 != "Excitatory neurons, midbrain"]
full_ann_sub = full_ann[full_ann["pbmc.merged.final@meta.data$annotations_1000"].str.contains("Excitatory neurons, midbrain")==False]
full_ann_sub

# full_ann = full_ann.set_index('id')
# glut_ann = glut_ann.set_index('id')


full_ann_final = pd.concat([full_ann_sub, glut_ann], axis=0, join='outer') 
full_ann_final

full_ann_t = full_ann.transpose()
glut_ann_t = glut_ann.transpose()


def return_indices_of_a(a, b):
  b_set = set(b)
  return [i for i, v in enumerate(a) if v in b_set]

################################################################################
p = "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/"
# iterate over rows in full_ann and then over glut_ann. 
# if rowname of glut_ann found in full_ann, get the index in full_ann.
# then replace the full_anns element with glut_ann's element
from collections import Counter

full_ann_final = full_ann.copy()
for cell_i, cell in enumerate(glut_ann.index):
    #print(cell_i)
    if cell in full_ann.index:
        i = list(full_ann.index).index(cell)
        # print(i)
        full_ann_final.iloc[i,0] = glut_ann.iloc[cell_i,0]
################################################################################
# full_ann_final.reset_index(inplace=True)
full_ann_final.rename(columns = {'pbmc.merged.final@meta.data$annotations_1000':'final_annotations'}, inplace = True)

for i, cell in enumerate(full_ann_final.iloc[:,0]):
    if cell == "Substantia nigra neurons":
        full_ann_final.iloc[i,0] = "Inhibitory neurons, midbrain"
        

Counter(list(full_ann_final.iloc[:,0]))

full_ann_final.to_csv(p + "final_annotations.csv")


full_ann_final = pd.merge(full_ann_t, glut_ann_t, how="outer")

full_ann_final = pd.merge(full_ann_t, glut_ann_t, left_index=True, right_index=True, how="inner")

