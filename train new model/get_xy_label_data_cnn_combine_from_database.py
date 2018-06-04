# Usage: python get_xy_data_cnn_combine_from_database.py bulk_gene_list.txt sc_gene_list.txt gene_pair_list  data_separation index list  bulk expression data  sc exprsssion data
# command line in developer's linux machine :
# python get_xy_label_data_cnn_combine_from_database.py bulk_gene_list.txt sc_gene_list.txt mmukegg_new_new_unique_rand_labelx.txt mmukegg_new_new_unique_rand_labelx_num.npy /home/yey3/sc_process_1/new_bulk_mouse/prs_calculation/mouse_bulk.h5 /home/yey3/sc_process_1/rank_total_gene_rpkm.h5
#################INPUT################################################################################################################################
# 1, bulk_gene_list.txt is the list that convert bulk expression data gene set into gene symbol IDs. format: 'gene symbol IDs\t bulk gene ID'
# 2, sc_gene_list.txt is the list that convert sc expression data gene set into gene symbol IDs. format: 'gene symbol IDs\t sc gene ID'
# 3, gene_pair_list is the list that contains gene pairs and their labels. format : 'GeneA    GeneB     0'
# 4, data_separation index list is a number list that divide gene_pair_list into small parts
# here we use data separation index list to divide gene pairs into small data parts, and make sure that the gene pairs in each index inteval is completely isolated from others. And we can evaluate CNNC's performance on only a small data part.
# if users do not need to separate data, they can just generate a index list to divide the data into N equal parts.
# 5, bulk expression data  it should be a hdf5 format. users can use their own data or data we provided.
# 6, sc expression data  it should be a hdf5 format. users can use their own data or data we provided.
#################OUTPUT
# it generate a data_label folder, and a series of data files containing Nxdata_tf (NEPDF file), ydata_tf (label file) and zdata_tf (gene symbol pair file) for each data part divided.

import pandas as pd
from numpy import *
import json, re,os, sys
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
        h[search_result.group(1).lower()]=search_result.group(2) # h [gene symbol] = gene ID
    s.close()
    return h


# Script starts from here
if len(sys.argv) < 7:
    print ('No enough input files')
    sys.exit()

h_gene_list_bulk =get_gene_list_bulk(sys.argv[1]) #'bulk_gene_list.txt')#
print ('read bulk gene list')
h_gene_list =get_gene_list(sys.argv[2]) # 'sc_gene_list.txt')#
print ('read sc gene list')



store = pd.HDFStore(sys.argv[5])#'/home/yey3/sc_process_1/new_bulk_mouse/prs_calculation/mouse_bulk.h5')  ### bulk RNA-seq expression data        )#
rpkm_bulk = store['rpkm']
store.close()
print ('read bulk RNA-seq expression')
store = pd.HDFStore(sys.argv[6])#'/home/yey3/sc_process_1/rank_total_gene_rpkm.h5')    # scRNA-seq expression data                        )#
rpkm = store['rpkm']
store.close()

########## generate NEPDF matrix
gene_pair_label = []
s=open(sys.argv[3])#'mmukegg_new_new_unique_rand_labelx.txt')#)   ### read the gene pair and label file
for line in s:
    gene_pair_label.append(line)
gene_pair_index = load(sys.argv[4])#'mmukegg_new_new_unique_rand_labelx_num.npy')#sys.argv[6]) # read file speration index
# ##### Here we saved data separation index as a np.array for convenience, users can also save it as a txt, and load it then convert it as a np.array
s.close()
gene_pair_label_array = array(gene_pair_label)
for i in range(len(gene_pair_index)-1):   #### many sperations
    print (i)
    start_index = gene_pair_index[i]
    end_index = gene_pair_index[i+1]
    x = []
    y = []
    z = []
    for gene_pair in gene_pair_label_array[start_index:end_index]: ## each speration
        separation = gene_pair.split()
        x_gene_name,y_gene_name,label = separation[0],separation[1],separation[2]
        x_tf_bulk = log10(rpkm_bulk[h_gene_list_bulk[x_gene_name]][0:249] + 10 ** -2)
        x_tf = log10(rpkm[int(h_gene_list[x_gene_name])][0:43261] + 10 ** -2)
        x_gene = log10(rpkm[int(h_gene_list[y_gene_name])][0:43261] + 10 ** -2)
        x_gene_bulk = log10(rpkm_bulk[h_gene_list_bulk[y_gene_name]][0:249] + 10 ** -2)
        z.append(x_gene_name+'\t'+y_gene_name)
        y.append(label)
        H_T = histogram2d(x_tf, x_gene, bins=32)
        H = H_T[0].T
        HT = (log10(H / 43261 + 10 ** -4) + 4) / 4
        H_T_bulk = histogram2d(x_tf_bulk, x_gene_bulk, bins=32)
        H_bulk= H_T_bulk[0].T
        HT_bulk = (log10(H_bulk / 43261 + 10 ** -4) + 4)/4
        x.append(concatenate((HT, HT_bulk), axis=0))
    if (len(x)>0):
        xx = array(x)[:, :, :, newaxis]
    else:
        xx = array(x)
    save(save_dir+'/Nxdata_tf' + str(i) + '.npy', xx)
    save(save_dir+'/ydata_tf' + str(i) + '.npy', array(y))
    save(save_dir+'/zdata_tf' + str(i) + '.npy', array(z))



