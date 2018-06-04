# CNNC
A supervised framework for co-expression analysis 
# Title: Convolutional neural network co-expression analysis (CNNC)

date: 2018-06-03

# tags:

# data sources
scRNA-seq :https://s3.amazonaws.com/scquery/processed_data/expr_data.hdf5

bulk RNA-seq : https://s3.us-east-2.amazonaws.com/mousebulkexprssion/mouse_bulk.h5

# — code

# code sources:

# Trained model for:

— mESC TF-target prediction model based on GTRD CHIP-seq database

— KEGG Pathway prediction model

— Reactome Pathway prediction model



# Train model for a new task:

# — readme

Using ‘trained models’, one can predict if one gene pair can interact as TF-target, KEGG pathway edges or Reactome protein interaction pair.

Using ‘train new model’,one can define a new predict task.

—- manual

# Trained model:

# step1, users need to provide gene pair candidate list;

# step2, use get_xy_data_cnn_combine_from_database.py to get gene pair NEPDF list;

#Usage: python get_xy_data_cnn_combine_from_database.py bulk_gene_list.txt sc_gene_list.txt gene_pair_list  data_separation_index list  bulk_expression_data  sc_exprsssion_data

#command line in developer's linux machine :

#python get_xy_data_cnn_combine_from_database.py bulk_gene_list.txt sc_gene_list.txt mmukegg_new_new_unique_rand_no_labelx.txt mmukegg_new_new_unique_rand_no_labelx_num.npy /home/yey3/sc_process_1/new_bulk_mouse/prs_calculation/mouse_bulk.h5 /home/yey3/sc_process_1/rank_total_gene_rpkm.h5

#################INPUT################################################################################################################################

#1, bulk_gene_list.txt is the list that convert bulk expression data gene set into gene symbol IDs. format: 'gene symbol IDs\t bulk gene ID'

#2, sc_gene_list.txt is the list that convert sc expression data gene set into gene symbol IDs. format: 'gene symbol IDs\t sc gene ID'

#3, gene_pair_list is the list that contains gene pairs and their labels. format : 'GeneA    GeneB   '

#4, data_separation index list is a number list that divide gene_pair_list into small parts

#here we use data separation index list to divide gene pairs into small data parts, and make sure that the gene pairs in each index inteval is completely isolated from others. And we can evaluate CNNC's performance on only a small data part.
#if users do not need to separate data, they can just generate a index list to divide the data into N equal parts.

#5, bulk expression data  it should be a hdf5 format. users can use their own data or data we provided.

#6, sc expression data  it should be a hdf5 format. users can use their own data or data we provided.

#################OUTPUT

#it generate a data_no_label folder, and a series of data files containing Nxdata_tf (NEPDF file)  and zdata_tf (gene symbol pair file) for each data part divided.

Here we use gene symbol information to align bulk, scRNA-seq and gene pair's gene sets. In our own data, scRNA-seq used entrez ID, bulk RNA-seq used ensembl ID, gene pair list used gene symbol ID, thus we used 'bulk_gene_list.txt' and 'sc_gene_list.txt' to convert all the IDs to gene symbols. Please make IDs convert to gene symbol ID files for bulk and scRNA-seq data.

# step3, use predict_no_y.py to do prediction;

#Usage: python predict_no_y.py number_of_data_parts_divided KEGG_or_GTRD_or_Reactome

#command line in developer's linux machine :

#python predict_no_y.py  2   KEGG

(In the models folder are three trained model for GTRD TF-target, KEGG and Reactome database respectively)

# Train new model:

# step1, users need to provide gene pair candidate list, their labels, and  and their expression data (optional);
gene pair candidate label list, such as mmukegg_new_new_unique_rand_labelx.txt
# step2, use get_xy_label_data_cnn_combine_from_database.py to get gene pair NEPDF list and their labels;

#Usage: python get_xy_data_cnn_combine_from_database.py bulk_gene_list.txt sc_gene_list.txt gene_pair_list  data_separation index list  bulk expression data  sc exprsssion data

#command line in developer's linux machine :

#python get_xy_label_data_cnn_combine_from_database.py bulk_gene_list.txt sc_gene_list.txt mmukegg_new_new_unique_rand_labelx.txt mmukegg_new_new_unique_rand_labelx_num.npy /home/yey3/sc_process_1/new_bulk_mouse/prs_calculation/mouse_bulk.h5 /home/yey3/sc_process_1/rank_total_gene_rpkm.h5

#################INPUT################################################################################################################################

#1, bulk_gene_list.txt is the list that convert bulk expression data gene set into gene symbol IDs. format: 'gene symbol IDs\t bulk gene ID'

#2, sc_gene_list.txt is the list that convert sc expression data gene set into gene symbol IDs. format: 'gene symbol IDs\t sc gene ID'

#3, gene_pair_list is the list that contains gene pairs and their labels. format : 'GeneA    GeneB     0'

#4, data_separation index list is a number list that divide gene_pair_list into small parts

#here we use data separation index list to divide gene pairs into small data parts, and make sure that the gene pairs in each index inteval is completely isolated from others. And we can evaluate CNNC's performance on only a small data part.
#if users do not need to separate data, they can just generate a index list to divide the data into N equal parts.

#5, bulk expression data  it should be a hdf5 format. users can use their own data or data we provided.

#6, sc expression data  it should be a hdf5 format. users can use their own data or data we provided.

#################OUTPUT

#it generate a data_label folder, and a series of data files containing Nxdata_tf (NEPDF file), ydata_tf (label file) and zdata_tf (gene symbol pair file) for each data part divided.


# step3, use train_with_labels_three_fold.py to train a new model with three-fold cross validation;

#Usage  python train_with_labels_wholedata.py number_of_data_parts_divided

#command line in developer's linux machine :

#module load cuda-8.0 using GPU

#srun -p gpu --gres=gpu:1 -c 2 --mem=20Gb python train_wtih_labels_wholedata.py 3057 > results_whole.txt

#######################OUTPUT

#it will generate a folder 'wholeXXXXX', in which 'keras_cnn_trained_model_shallow.h5' is the final trained model

# step4, use train_with_labels_wholedata.py to train a new model with whole data;

#Usage  python train_with_labels_wholedata.py number_of_data_parts_divided

#command line in developer's linux machine :

#module load cuda-8.0 using GPU

#srun -p gpu --gres=gpu:1 -c 2 --mem=20Gb python train_wtih_labels_wholedata.py 3057 > results_whole.txt

#######################OUTPUT

#it will generate a folder 'wholeXXXXX', in which 'keras_cnn_trained_model_shallow.h5' is the final trained model

