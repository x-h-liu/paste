from operator import index
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns

# gene_expression_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice1/filtered_feature_bc_matrix.h5"
# adata = sc.read_10x_h5(gene_expression_file)
# adata.var_names_make_unique()
# print(adata.shape)
# # print(adata.obs.index)

# location_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice1/spatial/tissue_positions_list.csv"
# all_locations = pd.read_csv(filepath_or_buffer=location_file, header=None, index_col=0, usecols=[0, 1, 4, 5])
# tissue_locations = all_locations[all_locations[1] == 1]
# spatial = []
# for spotname in adata.obs.index:
#     spatial.append([tissue_locations.loc[spotname][5], tissue_locations.loc[spotname][4]])
# adata.obsm['spatial'] = np.array(spatial)

# print(adata.shape)

def plot_slice(slice_,figsize=None,ax=None, s=100):
    (min_x,min_y),(max_x,max_y) = slice_.obsm['spatial'].min(axis=0),slice_.obsm['spatial'].max(axis=0)
    len_x,len_y=max_x-min_x,max_y-min_y
    if not figsize: figsize=(10*(len_x/max(len_x,len_y)),10*(len_y/max(len_x,len_y)))
    if not ax: plt.figure(figsize=figsize)
    g = sns.scatterplot(x=slice_.obsm['spatial'][:,0],y=slice_.obsm['spatial'][:,1],linewidth=0,s=s, marker=".",ax=ax)
    if not ax: ax=g
    if ax:
        ax.invert_yaxis()
        ax.axis('off')
    plt.show()


def plot_slice_umi(slice_, figsize=None, ax=None, s=100):
    (min_x,min_y),(max_x,max_y) = slice_.obsm['spatial'].min(axis=0), slice_.obsm['spatial'].max(axis=0)
    len_x, len_y = max_x-min_x, max_y-min_y
    if not figsize: figsize=(10*(len_x/max(len_x,len_y)),10*(len_y/max(len_x,len_y)))
    if not ax: plt.figure(figsize=figsize)
    min_umi = np.min(slice_.obs["sum_umi"])
    max_umi = np.max(slice_.obs["sum_umi"])
    norm = colors.Normalize(vmin=min_umi, vmax=max_umi, clip=True)
    # cmap = sns.color_palette("rocket_r", as_cmap=True)
    cmap = sns.color_palette("YlOrBr", as_cmap=True)
    c = [cmap(norm(i)) for i in slice_.obs["sum_umi"]]
    g = sns.scatterplot(x = slice_.obsm['spatial'][:,0],y = slice_.obsm['spatial'][:,1],linewidth=0,s=s, marker=".", c=c, ax=ax)
    if not ax: ax=g
    if ax:
        ax.invert_yaxis()
        ax.axis('off')


def read_HT(gene_expression_file, location_file, tumor_annotation_file):
    adata = sc.read_10x_h5(gene_expression_file)
    adata.var_names_make_unique()

    all_locations = pd.read_csv(filepath_or_buffer=location_file, header=None, index_col=0, usecols=[0, 1, 4, 5])
    tissue_locations = all_locations[all_locations[1] == 1]
    spatial = []
    for spotname in adata.obs.index:
        spatial.append([tissue_locations.loc[spotname][5], tissue_locations.loc[spotname][4]])
    adata.obsm['spatial'] = np.array(spatial)
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    gene_expression_matrix = adata.X
    spot_umi_counts = np.sum(gene_expression_matrix, axis=1)
    adata.obs['sum_umi'] = spot_umi_counts

    tumor_annotation_df = pd.read_csv(tumor_annotation_file, sep='\t', header=0, index_col=0)
    tumor_annotation = []
    for spotname in adata.obs.index:
        tumor_annotation.append(tumor_annotation_df.loc[spotname]['annotation'])
    adata.obs['tumor'] = tumor_annotation
    return adata


if __name__ == "__main__":
    # # HT 225 slice 1
    # # gene_expression_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice1/filtered_feature_bc_matrix.h5"
    # # location_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice1/spatial/tissue_positions_list.csv"

    # HT 225 slice 2
    # gene_expression_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice2/filtered_feature_bc_matrix.h5"
    # location_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice2/spatial/tissue_positions_list.csv"

    # HT 225 slice 3
    # gene_expression_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice3/filtered_feature_bc_matrix.h5"
    # location_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice3/spatial/tissue_positions_list.csv"

    # # HT 225 slice 4
    # # gene_expression_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice4/filtered_feature_bc_matrix.h5"
    # # location_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice4/spatial/tissue_positions_list.csv"

    # # HT 225 slice 5
    # # gene_expression_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice5/filtered_feature_bc_matrix.h5"
    # # location_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice5/spatial/tissue_positions_list.csv"

    # # HT 112
    # gene_expression_file = "/n/fs/ragr-data/datasets/DingLab/HT112/U2/filtered_feature_bc_matrix.h5"
    # location_file = "/n/fs/ragr-data/datasets/DingLab/HT112/U2/spatial/tissue_positions_list.csv"

    # adata = read_HT(gene_expression_file, location_file)
    # spot_umi_counts = adata.obs['sum_umi']
    # print(adata.shape)
    # print(spot_umi_counts.shape)
    # print("mean UMI per spot is: " + str(np.mean(spot_umi_counts)))
    # print("median UMI per spot is: " + str(np.median(spot_umi_counts)))
    # print("max UMI per spot is: " + str(np.max(spot_umi_counts)))
    # print("min UMI per spot is: " + str(np.min(spot_umi_counts)))

    # plot_slice(adata)
    # plot_slice_umi(adata)

    # counts, edges, bars = plt.hist(adata.obs['sum_umi'])
    # plt.bar_label(bars)
    # plt.xlabel("UMI")
    # plt.ylabel("Number of spots")
    # plt.show()


    tumor_divide_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice1/HT225C1-Th1_ST.metadata.tsv"
    tumor_divide = pd.read_csv(tumor_divide_file, sep='\t', header=0, index_col=0)
    print(tumor_divide)
    print(tumor_divide.loc['CAAGCAACGTCGGAGT-1']['annotation'])


