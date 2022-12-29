import os
import pandas as pd
import anndata as ad
import numpy as np
#make gene-cell binary matrix 
genes = ["TH", "DDC", "GCH1", "VMAT", "SLC6A3","EPHA4", "CHRNA5", "NRIP3", 
                                    "KCNS3", "CPLX1", "SOX6", "NDNF", "ALDH11A", "ALDH1A7",
                                    "LIX1", "PTPNS", "NOSTRIN", "SERPINE2", "KCNIP3",
          "FJX1", "FOXA2", "EN2", "NTF3", "GFRA2", "KLFC3", "CALB1", "CHST8", "CLBN4", "LPL", 
        "NHKH2", "OTX1", "OTX2", "NEUROD6", "GRP", "TCF12", "CALCA", "NTF3", "GFRA2", "VIP", "CCK",
         "SYN2", "CBLN1", "GPX3", "SLC32A1", "CTXN3", "ETC1", "LMX1A"]

genes = list(np.unique(np.asarray(genes)))

cell_types_da = {'DA_SNc':["TH", "DDC", "GCH1", "SLC18A2", "SLC6A3", "EPHA4", "CHRNA5", "NRIP3", "KCNS3", "CPLX1", "SOX6", "NDNF", "ALDH11A", "ALDH1A7",
"LIX1", "PTPNS", "NOSTRIN", "SERPINE2", "KCNIP3"], 'DA_VTA1':["TH", "DDC", "GCH1", "VMAT", "SLC6A3", "EPHA4", "CHRNA5", "NRIP3", "KCNS3", "CPLX1", "SOX6", "NDNF",
             "FJX1", "FOXA2", "EN2", "NTF3", "GFRA2"], 'DA_VTA2':["TH", "DDC", "GCH1", "VMAT", "SLC6A3", "KLFC3", "CALB1", "CHST8", "CLBN4", "LPL", "ALDH11A", "ALDH1A7", 
            "NHKH2", "OTX1", "OTX2", "NEUROD6", "GRP", "TCF12", "CALCA"],
            'DA_VTA3':["TH", "DDC", "GCH1", "VMAT", "SLC6A3", "KLFC3", "CALB1", "CHST8", "CBLN4", "LPL", "NHLH2", "OTX1", "SYN2",
             "CBLN1", "GPX3", "FJX1", "FOXA2", "EN2","NTF3", "GFRA2", "VIP", "CCK"],
             'DA_VTA4':["TH", "DDC", "GCH1", "VMAT", "SLC6A3", "KLFC3", "CALB1", "CHST8", "SYN2", "CBLN1","GPX3", "SLC32A1", "CTXN3", "ETC1", "LMX1A"]}

len(cell_types_da['DA_SNc'])
len(cell_types_da['DA_VTA1'])
len(cell_types_da['DA_VTA2'])
len(cell_types_da['DA_VTA3'])


cells.keys()
genes = ["TH", "DDC", "GCH1", "SLC18A2", "ALDH1A1", "SOX6", "VLGUT2", "OTX2", "VGAT", "VIP", "CCK", "CALB1" ]
cells = {'Aldh1a1+/Sox6+':["TH", "DDC", "GCH1", "SLC18A2", "ALDH1A1", "SOX6"], 
        'Aldh1a1+/Sox6-':["TH", "DDC", "GCH1", "SLC18A2", "ALDH1A1"], 
        'Vglut2+/Aldh1a1- (SNc)':["TH", "DDC", "GCH1", "SLC18A2", "VLGUT2"], 
        'Aldh1a1+/Otx2+/Vglut2+': ["TH", "DDC", "GCH1", "SLC18A2", "VLGUT2", "ALDH1A1","OTX2"],
        'Vglut2+/Aldh1a1- (VTA)': ["TH", "DDC", "GCH1", "SLC18A2", "VLGUT2"],
        'Vgat+/Vglut2': ["TH", "DDC", "GCH1", "SLC18A2", "VLGUT2", "VGAT"],
        'VIP/VGLUT2':["TH", "DDC", "GCH1", "SLC18A2", "VLGUT2", "VIP"],
        'Poorly defined': ["TH", "DDC", "GCH1", "SLC18A2", "VLGUT2", "CCK", "CALB1"]}

genes = ["TH", "DDC", "GCH1", "SLC18A2", "ALDH1A1", "SOX6", "OTX2", "VIP", "CCK", "CALB1" ]

cell_types_da = {'Aldh1a1+/Sox6+':["TH", "DDC", "GCH1", "SLC18A2", "ALDH1A1", "SOX6"], 
        'Aldh1a1+/Sox6-':["TH", "DDC", "GCH1", "SLC18A2", "ALDH1A1"], 
        'Vglut2+/Aldh1a1- (SNc)':["TH", "DDC", "GCH1", "SLC18A2"], 
        'Aldh1a1+/Otx2+/Vglut2+': ["TH", "DDC", "GCH1", "SLC18A2", "ALDH1A1","OTX2"],
        'Vglut2+/Aldh1a1- (VTA)': ["TH", "DDC", "GCH1", "SLC18A2"],
        'VIP/VGLUT2':["TH", "DDC", "GCH1", "SLC18A2", "VIP"],
        'Poorly defined': ["TH", "DDC", "GCH1", "SLC18A2", "CCK", "CALB1"],
        'DA_SNc':["TH", "DDC", "GCH1", "SLC18A2", "SLC6A3", "EPHA4", "CHRNA5", "NRIP3", "KCNS3", "CPLX1", "SOX6", "NDNF", "ALDH11A", "ALDH1A7",
"LIX1", "PTPNS", "NOSTRIN", "SERPINE2", "KCNIP3"], 'DA_VTA1':["TH", "DDC", "GCH1", "VMAT", "SLC6A3", "EPHA4", "CHRNA5", "NRIP3", "KCNS3", "CPLX1", "SOX6", "NDNF",
             "FJX1", "FOXA2", "EN2", "NTF3", "GFRA2"], 'DA_VTA2':["TH", "DDC", "GCH1", "VMAT", "SLC6A3", "KLFC3", "CALB1", "CHST8", "CLBN4", "LPL", "ALDH11A", "ALDH1A7", 
            "NHKH2", "OTX1", "OTX2", "NEUROD6", "GRP", "TCF12", "CALCA"],
            'DA_VTA3':["TH", "DDC", "GCH1", "VMAT", "SLC6A3", "KLFC3", "CALB1", "CHST8", "CBLN4", "LPL", "NHLH2", "OTX1", "SYN2",
             "CBLN1", "GPX3", "FJX1", "FOXA2", "EN2","NTF3", "GFRA2", "VIP", "CCK"],
             'DA_VTA4':["TH", "DDC", "GCH1", "VMAT", "SLC6A3", "KLFC3", "CALB1", "CHST8", "SYN2", "CBLN1","GPX3", "SLC32A1", "CTXN3", "ETC1", "LMX1A"]}
genes = ["TH", "DDC", "GCH1", "VMAT", "SLC6A3","EPHA4", "CHRNA5", "NRIP3", 
                                    "KCNS3", "CPLX1", "SOX6", "NDNF", "ALDH11A", "ALDH1A7",
                                    "LIX1", "PTPNS", "NOSTRIN", "SERPINE2", "KCNIP3",
          "FJX1", "FOXA2", "EN2", "NTF3", "GFRA2", "KLFC3", "CALB1", "CHST8", "CLBN4", "LPL", 
        "NHKH2", "OTX1", "OTX2", "NEUROD6", "GRP", "TCF12", "CALCA", "NTF3", "GFRA2", "VIP", "CCK",
         "SYN2", "CBLN1", "GPX3", "SLC32A1", "CTXN3", "ETC1", "LMX1A"]

def reformat(name):
    return(name.upper())

ref = pd.read_csv('/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/rds_objects/Integ_count_matrix.csv')
ref = pd.read_csv("/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/da_subtype_counts.csv")
ref = pd.read_csv("/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/annotations_final_nonspecific_da_only.csv")

ref = pd.read_csv("/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/da_subtype_counts_scifi.csv")


ref.rename(columns = {'Unnamed: 0':'genes'}, inplace = True)
ref.index = ref['genes']
ref = ref.drop('genes', axis=1)
gene_row_upper = list(map(reformat, list(ref.index)))
ref.index = gene_row_upper

#genes = pd.read_csv('/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/rds_objects/Integ_genes.csv')
# ref




def determine_cell_type(cell_type_genes, df_da, set_th: float):
    match_info = {}  #key : da cell type, value: [index, perc of cells matching the ref genes]. for each cell, you need to check if
    #the cell, i.e the cell index has already been labelled prior, if it has, you need to see if the current one has a higher perc 
    #match compared to the prev one
    da_types_cell={}
    undefined = []
    for da_type in cell_type_genes.keys():
        da_type_coord = []
        match_data = []
        for col_i, cell in enumerate(df_da.columns):
            found_genes_cell=[]
            for row_i in range(len(df_da.index)-1):
                if df_da.index[row_i] in cell_type_genes[da_type]:
                    if float(df_da.iloc[row_i,col_i])>0:
                        found_genes_cell.append(df_da.index[row_i])
                        # print(f"match {df_da.index[row_i]} found")
            #print(found_genes_cell)
            # check that TH is found in the list as it should exist in all DA cells
            if "TH" and "SLC6A3" in found_genes_cell:
                if len(found_genes_cell)/len(cell_type_genes[da_type])>=set_th:
                    perc = len(found_genes_cell)/len(cell_type_genes[da_type])
                    print(f'Found genes per cell: {len(found_genes_cell)} with perc: {perc}')
                    da_type_coord.append(int(col_i))
                    da_type_coord.append(perc)
                    # round(len(found_genes_cell)/len(cell_type_genes), 6)
                    # match_data.append([col_i, round(len(found_genes_cell)/len(cell_type_genes), 6)])
                    # print([col_i, round(len(found_genes_cell)/len(cell_type_genes), 6)])
                #da_types_cell['DA_A9'].append(col_i)
        # match_info[da_type]=match_data
            else:
                # undefined cells can be found in multiple cell types. thus, add the list into the key dict
                # at the very end, not after iterating through an individual da type
                print("Unable to define the cell type as DA, setting as Undefined...")
                perc = len(found_genes_cell)/len(cell_type_genes[da_type])
                da_type_coord.append(int(col_i))
                da_type_coord.append(perc)
                undefined.append(da_type_coord)
        da_types_cell[da_type] = da_type_coord
    da_types_cell['Neurons undefined'] = da_type_coord

    #da_types_cell[cell_type]=da_type_coord
    #da_dict_ann = clean_da_annot(match_info, da_types_cell)
    return(da_types_cell)
ref.columns

data_meta = determine_cell_type(cell_types_da, ref, 0.5)
data_meta1 = determine_cell_type(cell_types_da, ref, 0.5)
data_meta

len(data_meta['Neurons undefined'])
len(data_meta['DA_SNc'])
len(data_meta['DA_VTA3'])


data_meta.keys()
data_meta['DA_SNc'][::2][-5:]
data_meta['DA_VTA1'][::2][-5:]


def da_subtype_annotation(ref, af):
    if 'Unnamed: 0' in ref.columns:
        ref.rename(columns = {'Unnamed: 0':'genes'}, inplace = True)
        ref.index = ref['genes']
        ref = ref.drop('genes', axis=1)
    gene_row_upper = list(map(reformat, list(ref.index)))
    ref.index = gene_row_upper
    data_meta = determine_cell_type(cell_types_da, ref, 0.5)
    data_meta1 = determine_cell_type(cell_types_da, ref, 0.5)
    da_type_dict = {}
    for da1 in data_meta1.keys():
        for da2 in data_meta1.keys():
            if da1 != da2:
                for perc_i, cell_i in zip(data_meta1[da1][1::2], data_meta1[da1][0::2]):
                    # cell_id_multiples = []
                    cell_id_multiples = {}
                    print(perc_i)
                    if cell_i in data_meta1[da2][::2]:
                        # cell_id_multiples.append([da1, perc_i])
                        cell_id_multiples[da1] = perc_i
                        # get the match index from data[da2] and get the perc based on that (the next element)
                        # cell_id_multiples.append([da2, data[da2][data[da2].index(cell_i) + 1]])
                        cell_id_multiples[da2] = data_meta1[da2][data_meta1[da2].index(cell_i) + 1]
                        if len(cell_id_multiples) > 0:
                            da_type_dict[cell_i] = cell_id_multiples

    for cell_duplicate in da_type_dict.keys():
        percs_cell_id = []
        print(f"===={cell_duplicate}=====")
        losing_da_type = min(da_type_dict[cell_duplicate], key=da_type_dict[cell_duplicate].get)
        winning_da_type = max(da_type_dict[cell_duplicate], key=da_type_dict[cell_duplicate].get)
        print(f'losing da type is {losing_da_type} with a perc of {da_type_dict[cell_duplicate][losing_da_type]} removing the cell id from the meta data')
        print("Cell duplicate: " + str(data_meta1[losing_da_type].index(cell_duplicate)))

        # if the difference between the winning and losing cell type is less than the ambiguity factor af, then
        # assign as a general DA cell
        # if da_type_dict[cell_duplicate][winning_da_type] - da_type_dict[cell_duplicate][losing_da_type] < af:
        #     data_meta1['Dopaminergic neurons, ventral midbrain (SNc, VTA)'] = data_meta1[losing_da_type][data_meta1[losing_da_type].index(cell_duplicate) + 1]
        cell_dups_perc = data_meta1[losing_da_type][data_meta1[losing_da_type].index(cell_duplicate) + 1]
        print(cell_dups_perc)
        # also remove the perc value next to the cell id
        data_meta1[losing_da_type].remove(cell_duplicate)
        data_meta1[losing_da_type].remove(cell_dups_perc)


    da_subtypes = pd.DataFrame(0, index=[0], columns=ref.columns).transpose()

    for da in data_meta1.keys():
        print(da)
        for cell_i, cell in enumerate(da_subtypes.index):
            print(cell_i)
            if cell_i in data_meta1[da]:
                da_subtypes.iloc[cell_i,0] = da       

    for i, cell_ident in enumerate(da_subtypes.iloc[:,0]):
        if cell_ident == 0:
            da_subtypes.iloc[i,0] = "Dopaminergic neurons, ventral midbrain (SNc, VTA)"
    return(da_subtypes)

da_subtypes = da_subtype_annotation(ref)
da_subtypes.to_csv("/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/DA_subtypes_parse_meta.csv")


'''
Find duplicates of index values via list comprehensions
Take up the given index value as the tuple (index, percentage)
compare between the samples
the one with the highest percentage will gain the identity
'''
data.keys()
'''
get the cell indeces, check if exists in the other da subtype. if it does, save into list of sub dicts
where the key is the cell id and the values are tuples containing the cell type + the given percentage
after all cell types have been gone through, go over the dict and see which cell type - perc pair
has the largest value. let the original dict keep this cell's index in its da type while the rest 
will remove theirs, including the value next to each cell id which is the perc.
once this list is ready, i
'''
data_meta['DA_VTA1']

data_meta1 = data_meta



da_type_dict = {}
for da1 in data_meta1.keys():
    for da2 in data_meta1.keys():
        if da1 != da2:
            for perc_i, cell_i in zip(data_meta1[da1][1::2], data_meta1[da1][0::2]):
                # cell_id_multiples = []
                cell_id_multiples = {}
                print(perc_i)
                if cell_i in data_meta1[da2][::2]:
                    # cell_id_multiples.append([da1, perc_i])
                    cell_id_multiples[da1] = perc_i
                    # get the match index from data[da2] and get the perc based on that (the next element)
                    # cell_id_multiples.append([da2, data[da2][data[da2].index(cell_i) + 1]])
                    cell_id_multiples[da2] = data_meta1[da2][data_meta1[da2].index(cell_i) + 1]
                    if len(cell_id_multiples) > 0:
                        da_type_dict[cell_i] = cell_id_multiples


'''
Go over each cell duplicate. each da_type_dict key is such a cell which contains the subdict with keys
of cell types that contain this cell id duplicate and the value being the percentage. get the da type 
which has the highest perc
the losing cell types get the general name of Dopaminergic neurons, ventral midbrain (SNc, VTA)

Update: subtract the highest vs lowest pairs. if the delta is more than a given threshold, define the DA
as the Dopaminergic neurons, ventral midbrain (SNc, VTA)
'''

len(da_type_dict)

winning_da_type = max(da_type_dict[1279], key=da_type_dict[1279].get)

i = data_meta1['DA_SNc'].index(7)

a = [i for i in data_meta1['DA_SNc'] if i ==7]
a

# data_meta1['DA_SNc'][data_meta1['DA_SNc'].index(7)]
data_meta1 = data_meta

def remove_datapoint_da(data_meta, cell_duplicate, da_type):
    try:
        cell_dups_perc = data_meta[da_type][data_meta[losing_da_type].index(cell_duplicate) + 1]
        # also remove the perc value next to the cell id
        data_meta[da_type].remove(cell_duplicate)
        data_meta[da_type].remove(cell_dups_perc)
        return(data_meta)
    except IndexError:
        return(data_meta)



da_type_dict[489]


# some duplicates not found in the original meta data. This may be because they got removed prior already

general_da_pool = []
for cell_duplicate in da_type_dict.keys():
    percs_cell_id = []
    print(f"===={cell_duplicate}=====")
    losing_da_type = min(da_type_dict[cell_duplicate], key=da_type_dict[cell_duplicate].get)
    winning_da_type = max(da_type_dict[cell_duplicate], key=da_type_dict[cell_duplicate].get)

    # if the difference between the winning and losing cell type is less than the ambiguity factor af, then
    # assign as a general DA cell. i.e. remove from the standard annotations (SNc, VTA...)
    if da_type_dict[cell_duplicate][winning_da_type] - da_type_dict[cell_duplicate][losing_da_type] < 0.1:
        print("Removing data point and assigning into the original pool of DA types")
        try:
            general_da_pool.append(data_meta1[winning_da_type][data_meta1[winning_da_type].index(cell_duplicate) + 1])
            data_meta1 = remove_datapoint_da(data_meta, cell_duplicate, winning_da_type)
        except ValueError:
            print(f"Duplicate {cell_duplicate} has already been removed from the cell type {winning_da_type}")
        try:
            general_da_pool.append(data_meta1[losing_da_type][data_meta1[losing_da_type].index(cell_duplicate) + 1])
            data_meta1 = remove_datapoint_da(data_meta, cell_duplicate, losing_da_type)
        except ValueError:
            print(f"Duplicate {cell_duplicate} has already been removed from the cell type {losing_da_type}")

    else:
        print(f'losing da type is {losing_da_type} with a perc of {da_type_dict[cell_duplicate][losing_da_type]} removing the cell id from the meta data')
        try:
            data_meta1 = remove_datapoint_da(data_meta, cell_duplicate, losing_da_type)
        except ValueError:
            print(f"Duplicate {cell_duplicate} has already been removed from the cell type {losing_da_type}")
data_meta1['Dopaminergic neurons, ventral midbrain (SNc, VTA)'] = general_da_pool
len(data_meta1['Dopaminergic neurons, ventral midbrain (SNc, VTA)'])

for cell_duplicate in da_type_dict.keys():
    percs_cell_id = []
    print(f"===={cell_duplicate}=====")
    losing_da_type = min(da_type_dict[cell_duplicate], key=da_type_dict[cell_duplicate].get)
    print(f'losing da type is {losing_da_type} with a perc of {da_type_dict[cell_duplicate][losing_da_type]} removing the cell id from the meta data')
    print("Cell duplicate: " + str(data_meta1[losing_da_type].index(cell_duplicate)))
    cell_dups_perc = data_meta1[losing_da_type][data_meta1[losing_da_type].index(cell_duplicate) + 1]
    print(cell_dups_perc)
    # also remove the perc value next to the cell id
    data_meta1[losing_da_type].remove(cell_duplicate)
    data_meta1[losing_da_type].remove(cell_dups_perc)

data_meta['DA_SNc'][::2][-5:]
data_meta['DA_VTA1'][::2][-5:]
len(data_meta['DA_SNc'])
len(data_meta1['DA_SNc'])

data_meta1["DA_VTA3"]
data_meta["DA_VTA3"]

data_meta["DA_SNc"]

data

da_subtypes = pd.DataFrame(0, index=[0], columns=ref.columns).transpose()

for da in data_meta1.keys():
    print(da)
    for cell_i, cell in enumerate(da_subtypes.index):
        print(cell_i)
        if cell_i in data_meta1[da]:
            da_subtypes.iloc[cell_i,0] = da       

da_subtypes
for i, cell_ident in enumerate(da_subtypes.iloc[:,0]):
    if cell_ident == 0:
        # da_subtypes.iloc[i,0] = "Dopaminergic neurons, ventral midbrain (SNc, VTA)"
        da_subtypes.iloc[i,0] = "Neurons undefined"

# da_subtypes_meta = da_subtypes.transpose()
da_subtypes.to_csv("/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/DA_subtypes_scifi_meta_2.csv")
# da_subtypes.to_csv("/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/DA_subtypes_parse_meta_2_harsh.csv")

l = {'a': [1,2,3]}
l['a'][2]

'''
now we need to remove the losing cell type from the original dict called data.
'''
data['DA_SNc'][]
res = [tup[0] for tup[0] in zip(data[1]['DA_SNc']) if tup[0] in data[0]['DA_VTA1']]

for a, b in zip(data[1]['DA_SNc']):
    print(a)



for tup in data[1]['DA_SNc']:
    print(tup)


data[0].keys()
data[1]
data[0]['DA_SNc']
data[0]['DA_VTA1']
# total = len(data[0]['Aldh1a1+/Sox6+']) + len(data[0]['Aldh1a1+/Sox6-']) + len(data[0]['Vglut2+/Aldh1a1- (SNc)']) + len(data[0]['Aldh1a1+/Otx2+/Vglut2+']) \
#     + len(data[0]['Vglut2+/Aldh1a1- (VTA)']) + len(data[0]['VIP/VGLUT2']) + len(data[0]['DA_SNc'])
total = len(data[0]['DA_SNc']) + len(data[0]['DA_VTA1']) + len(data[0]['DA_VTA2']) + len(data[0]['DA_VTA3']) \
    + len(data[0]['DA_VTA4'])
total

# da_subtypes = da_subtypes.rename_axis('cells').reset_index()

if 1691 in data[0]['DA_VTA1']:
    print(True)

da_subtypes.reset_index(inplace=True)
da_subtypes = ref.rename(columns = {'index':'genes'})

da_subtypes = ref.iloc[:,2:]
da_subtypes['DA_subtype'] = 0

ref.reset_index(inplace=True)
ref = ref.rename(columns = {'index':'genes'})
ref.reset_index(inplace=True)
ref = ref.rename(columns = {'index':'gene_index'})
ref
#if one cell passes the threshold 

#da_types_cell={'DA_A9': [], 'DA_SNc' : [], 'DA_VTA1' : [], 'DA_VTA2' : [], 'DA_VTA3' : [], 'DA_VTA4' : []}
da_base = determine_cell_type(base_markers_da.values(), df, 0.1)
#da_types_cell = list(map(determine_cell_type, cell_types_da, df_da))

#get the indeces of the ones which have sufficient number of base DA genes
da_base = {k: determine_cell_type(v, df, 0.2) for k, v in base_markers_da.items()}
len(da_base['DA_A9'])


#also save info on the number of correct genes found per cell type (when iterating through a cell)
def determine_cell_type(cell_type_genes, df_da, set_th: float):
    da_type_coord = []
    match_info = {}  #key : da cell type, value: [index, perc of cells matching the ref genes]. for each cell, you need to check if
    #the cell, i.e the cell index has already been labelled prior, if it has, you need to see if the current one has a higher perc 
    #match compared to the prev one
    da_types_cell={}
    for col_i, cell in enumerate(df_da.columns):
        found_genes_cell=[]
        for row_i in range(len(df_da.index)-1):
            if df_da.index[row_i] in cell_type_genes:
                if float(df_da.iloc[row_i,col_i])>0:
                    found_genes_cell.append(df_da.index[row_i])
                    print(f"match {df_da.index[row_i]} found")
        #print(found_genes_cell)
        if len(found_genes_cell)/len(cell_type_genes)>=set_th:
            da_type_coord.append(int(col_i))
            match_info[]=[col_i, round(len(found_genes_cell)/len(cell_type_genes), 5)]
            #da_types_cell['DA_A9'].append(col_i)
    #da_types_cell[cell_type]=da_type_coord
    return(da_type_coord)

gene_cell_bin={}
gene_cell_scina={}

gene_cell_bin["Genes"]=genes
gene_per_cell = []
for cell_type, cell_gene in cells.items():
    print(cell_type)
    binary_count=[]
    for gene in genes:
        if gene in cell_gene:
            binary_count.append(1)
            gene_per_cell.append(gene)
        else:
            binary_count.append(0)
    gene_cell_bin[cell_type]=binary_count
    gene_cell_scina[cell_type]=gene_per_cell
gene_cell_bin.keys()

gene_cell_bin_df = pd.DataFrame.from_dict(gene_cell_scina)
gene_cell_bin_df.set_index("Genes")

gene_cell_bin_df

gene_cell_bin_df.to_csv(loom_dir+"/DA_ref.csv", index=False)


