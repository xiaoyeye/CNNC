import pandas as pd
from collections import defaultdict
import matplotlib, operator
def load_json (file):
    import json
    with open (file, "r") as f:
        hx = json.load(f)
    return hx

cell_types = load_json('/home/yey3/nn_project2/data/mesc_cnnc/specific_experiment_term_mapping.json')#experiment annotation
store = pd.HDFStore('/home/yey3/sc_process_1/true_expression.h5') ##sc RNA-seq data
x = store['rpkm']
store.close()
sample_name = x.index

h_cell_type = defaultdict(list)
h_cell_exp = defaultdict(list)
indexx = 0
for sample in sample_name:
    for cell_name in cell_types[sample]:
        xx = cell_name.split(' ')
        yy = sample.split('_')
        cell_namex = " ".join(xx[1:])
        h_cell_type[cell_namex].append(indexx)
        if yy[0] not in h_cell_exp[cell_namex]:
            h_cell_exp[cell_namex].append(yy[0])
    indexx = indexx+1

h_cell_type_num = {}
for i in h_cell_type.keys():
    h_cell_type_num[i] = len(h_cell_type[i])

sort_h_cell_type = sorted(h_cell_type_num.items(), key=operator.itemgetter(1),reverse=True)
s= open ('/home/yey3/nn_project2/data/mesc_cnnc/cell_sample_num.txt','w')
for i in sort_h_cell_type:
    s.write(i[0]+'\t'+str(i[1])+'\n')

s.close()

cell_type_list = ['bone marrow cell','dendritic cell','embryonic stem cell'] ###  specific cell types
cell_type_listx = ['bone_marrow_cell','dendritic_cell','embryonic_stem_cell']
for i in range(len(cell_type_list)):
    store = pd.HDFStore('/home/yey3/nn_project2/data/mesc_cnnc/'+cell_type_listx[i]+'.h5')
    store['RPKMs'] = x.iloc[h_cell_type[cell_type_list[i]],:]
    store.close()

