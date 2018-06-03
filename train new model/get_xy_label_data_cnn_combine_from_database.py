# Usage: python get_xy_data_cnn_combine_from_database.py bulk_gene_list.txt,sc_gene_list.txt,gene_pair list, bulk expression data , sc exprsssion data
# python get_xy_label_data_cnn_combine_from_database.py bulk_gene_list.txt sc_gene_list.txt mmukegg_new_new_unique_label_head.txt /home/yey3/sc_process_1/new_bulk_mouse/prs_calculation/mouse_bulk.h5 /home/yey3/sc_process_1/rank_total_gene_rpkm.h5
import pandas as pd
from numpy import *
import json
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# import seaborn as sns
import scipy.stats as stats
import re
import os,sys
save_dir = os.path.join(os.getcwd(),'data_label')
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)
def get_gene_list_bulk(file_name):
    import re
    h={}
    s = open(file_name,'r')   #gene symbol ID list of bulk RNA-seq
    for line in s:
        search_result = re.search(r'^([^\s]+)\s+([^\s]+)',line)
        h[search_result.group(1).lower()]=search_result.group(2)   # h [gene symbol] = gene ID
    s.close()
    return h

def get_gene_list(file_name):
    import re
    h={}
    s = open(file_name,'r') #gene symbol ID list of sc RNA-seq
    for line in s:
        search_result = re.search(r'^([^\s]+)\s+([^\s]+)',line)
        h[search_result.group(1).lower()]=search_result.group(2).lower() # h [gene symbol] = gene ID
    s.close()
    return h

def get_tf_target_chipseq(file_name):
    import re
    from collections import defaultdict
    h_tf_target = defaultdict(list)
    h_tf_target_label = {}
    s_tf_target = open(file_name, 'r')  ## gene pair list
    for line in s_tf_target:
        search_result = re.search(r'^([^\s]+)\s+([^\s]+)\s+(\d+)', line)
        if search_result:
            if search_result.group(3) == '1':
                h_tf_target[search_result.group(1).lower()].append(search_result.group(2).lower())
            h_tf_target_label[search_result.group(1).lower()+'\t'+search_result.group(2).lower()] = search_result.group(3)
    s_tf_target.close()
    return h_tf_target,h_tf_target_label

##########

# Script starts from here
if len(sys.argv) < 6:
    print ('No enough input files')
    sys.exit()

h_gene_list_bulk =get_gene_list_bulk(sys.argv[1])
print ('read bulk gene list')
h_gene_list =get_gene_list(sys.argv[2])
print ('read sc gene list')
h_tf_target,h_tf_target_label=get_tf_target_chipseq(sys.argv[3])
print ('read gene pair list')

store = pd.HDFStore(sys.argv[4])#'/home/yey3/sc_process_1/new_bulk_mouse/prs_calculation/mouse_bulk.h5')  ### bulk RNA-seq expression data
rpkm_bulk = store['rpkm']
store.close()
print ('read bulk RNA-seq expression')
store = pd.HDFStore(sys.argv[5])#'/home/yey3/sc_process_1/rank_total_gene_rpkm.h5')    # scRNA-seq expression data
rpkm = store['rpkm']
store.close()
print ('read scRNA-seq expression')
gene_list = list(h_gene_list.keys())
gene_list_bulk = list(h_gene_list_bulk.keys())
kegg_tf_list = list(h_tf_target.keys())
ovlp_tf_set = list(set(gene_list ) & set(gene_list_bulk) & set (kegg_tf_list ))
h_ovlp_tf_set = {}
h_ovlp_tf_set['ovlp'] = ovlp_tf_set
with open("ovlp_tf_set_new_unique.json",'w') as handle:
   json.dump(h_ovlp_tf_set,handle)

print ('total genes in position 0:',len(ovlp_tf_set))#the data was divided by genes in position 0.
#with open("/home/yey3/cnn_project/chip-seq/keggp/new/ovlp_tf_set_new_unique.json") as handle:
#    h_ovlp_tf_set = json.load(handle)
#ovlp_tf_set = h_ovlp_tf_set['ovlp']

########## generate NEPDF matrix
for i in range(len(ovlp_tf_set)):
    x = []
    y = []
    z = []
    x_tf_bulk = log10(rpkm_bulk[h_gene_list_bulk[ovlp_tf_set[i]]][0:249]+10**-2)
    x_tf = log10(rpkm[int(h_gene_list[ovlp_tf_set[i]])][0:43261]+10**-2)
    print (i)
    for j in h_tf_target[ovlp_tf_set[i]]:#(total):
        if j in h_gene_list and j in h_gene_list_bulk:
            x_gene = log10(rpkm[int(h_gene_list[j])][0:43261]+10**-2)
            x_gene_bulk = log10(rpkm_bulk[h_gene_list_bulk[j]][0:249]+10**-2)
            #y.append(1)  ### x regulate y
            #y.append(2)  ####x does not regulate y, instead y regulate x
            z.append(ovlp_tf_set[i]+'\t'+j)
            y.append(h_tf_target_label[ovlp_tf_set[i]+'\t'+j])
            #z.append(j+'\t'+ovlp_tf_set[i])
            H_T = histogram2d(x_tf, x_gene, bins=32)
            H= H_T[0].T
            HT = (log10(H / 43261 + 10 ** -4) + 4)/4
            HTT = HT.T
            H_T_bulk = histogram2d(x_tf_bulk, x_gene_bulk, bins=32)
            H_bulk= H_T_bulk[0].T
            HT_bulk = (log10(H_bulk / 43261 + 10 ** -4) + 4)/4
            HTT_bulk = HT_bulk.T
            x.append(concatenate((HT, HT_bulk), axis=0))
            #x.append(concatenate((HTT, HTT_bulk), axis=0))
    if (len(x)>0):
        xx = array(x)[:, :, :, newaxis]
    else:
        xx = array(x)
    save(save_dir+'/Nxdata_tf' + str(i) + '.npy', xx)
    save(save_dir+'/ydata_tf' + str(i) + '.npy', array(y))
    save(save_dir+'/zdata_tf' + str(i) + '.npy', array(z))



#plt.imshow(x[1], interpolation='nearest', origin='low')
#plt.savefig('/home/yey3/cnn_project/chip-seq/keggp/new/1.png')


