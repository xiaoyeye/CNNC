# CNNC
A supervised framework for co-expression analysis 
#Title: Convolutional neural network co-expression analysis (CNNC)

date: 2018-06-03

#tags:

#data sources
scRNA-seq :https://s3.amazonaws.com/scquery/processed_data/expr_data.hdf5

bulk RNA-seq : https://s3.us-east-2.amazonaws.com/mousebulkexprssion/mouse_bulk.h5

#— code

#code sources:

#Trained model for:
— mESC TF-target prediction model based on GTRD CHIP-seq database

— KEGG Pathway prediction model

— Reactome Pathway prediction model



#train model for a new task:

— readme

Using ‘Trained models’, one can predict if one gene pair can interact as TF-target, KEGG pathway edges or Reactome protein interaction pair.

Using ‘train new model’,one can define a new predict task.

—- manual

#Trained model:

#step1, users need to provide gene pair candidate list or expression data;

#step2, use get_xy_data_cnn_combine_from_database.py to get gene pair NEPDF list;

Usage: python get_xy_data_cnn_combine_from_database.py bulk_gene_list.txt, sc_gene_list.txt, gene_pair list, bulk expression data, sc exprsssion data

#step3, use predict_no_y.py to do prediction;

Usage: predict_no_y.py number of genes in position0, kEGG or GTRD or Reactome

(In the models are three trained model for GTRD TF-target, KEGG and Reactome database respectively)

#Train new model:

#step1, users need to provide gene pair candidate list and their labels;

#step2, use get_xy_label_data_cnn_combine_from_database.py to get gene pair NEPDF list and their labels;

Usage: python get_xy_data_cnn_combine_from_database.py bulk_gene_list.txt, sc_gene_list.txt, gene_pair list, bulk expression data, sc exprsssion data

#step3, use train_with_labels.py to train a new model;

Usage: python train_with_labels.py
